use crate::stats::{GlobalPi, GlobalStatistic, SiteComposable};
use crate::PopgenError;
use proptest::collection::vec;
use proptest::proptest;
use rand::distr::uniform::SampleRange;

proptest!(
// we don't need to invoke actual parallelism, just show that the component-wise approach is correct
#[test]
fn pi_equivalent_concurrent(
    seed in 0..u64::MAX,
    ploidy in 1_usize..5,
    num_samples in 2_usize..50,
    missing_data_rate_raw in 0_f64..1.0 - f64::EPSILON,
    num_sites in 1_usize..10,
    non_normalized_freqs in vec(vec(f64::EPSILON..1_f64, 1..10), 10),
) {
    use rand::prelude::*;

    let mut rng = StdRng::seed_from_u64(seed);
    let mut sites = vec![];
    for nnf in non_normalized_freqs.iter().take(num_sites) {
        let sum = nnf.iter().sum::<f64>();
        let freqs = nnf.iter().map(|fi| fi / sum).collect::<Vec<_>>();
        let missing_data_rate = if missing_data_rate_raw > 0.0 {
            Some(crate::testing::testdata::RandomSiteOptions {
                missing_data_rate: Some(missing_data_rate_raw),
            })
        } else {
            None
        };
        let site = crate::testing::testdata::random_site_rng(
            num_samples,
            ploidy,
            &freqs,
            missing_data_rate,
            &mut rng,
        );
        if !site.iter().flat_map(|i| i.iter()).all(|i| i.is_none()) {
            sites.push(site);
        }
    }
    let counts = crate::testing::testdata::single_pop_counts(&mut sites.iter());
    let pi = counts.iter().try_fold(GlobalPi::default(), |mut pi, s| {
        pi.try_add_site(s)?;
        Ok::<_, PopgenError>(pi)
    });

    let mut components = counts.iter().map(|s| {
        let mut pi = GlobalPi::default();
        pi.try_add_site(s)?;
        Ok::<_, PopgenError>(pi)
    }).collect::<Vec<_>>();

    // let's combine the components out of order and with random associativity
    // this test is to make sure that this is equivalent to the serial case
    components.shuffle(&mut rng);
    while components.len() > 1 {
        let pick1 = components.remove((..components.len()).sample_single(&mut rng).unwrap());
        let pick2 = components.remove((..components.len()).sample_single(&mut rng).unwrap());
        match (pick1, pick2) {
            (Err(one), Ok(_)) | (Ok(_), Err(one)) | (Err(one), Err(_)) => {
                components.push(Err(one));
            },
            (Ok(mut one), Ok(two)) => {
                components.push(one.try_combine(&two).map(|_| one));
            },
        }
    }

    let parallel = if !components.is_empty() {
        components.remove(0)
    } else {
        Ok(GlobalPi::default())
    };

    match pi {
        Err(e) => {
            assert_eq!(std::mem::discriminant(&e), std::mem::discriminant(&parallel.err().unwrap()));
        },
        Ok(pi) => {
            let pi_raw = pi.as_raw();
            let parallel_raw = parallel.unwrap().as_raw();
            assert!((pi_raw - parallel_raw).abs() < 0.00001)
        },
    }
}
);
