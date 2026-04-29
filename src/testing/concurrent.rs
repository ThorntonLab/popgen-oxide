use crate::stats::{GlobalPi, GlobalStatistic, SiteComposable, WattersonTheta};
use crate::testing::testdata::Site;
use crate::PopgenError;
use proptest::collection::vec;
use proptest::proptest;
use rand::distr::uniform::SampleRange;
use rand::rngs::StdRng;
use rand::seq::SliceRandom;
use rand::SeedableRng;

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
    let mut rng = StdRng::seed_from_u64(seed);
    let sites = super::testdata::make_random_sites(
        &mut rng,
        ploidy,
        num_samples,
        missing_data_rate_raw,
        num_sites,
        non_normalized_freqs,
    );

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

proptest!(
// we don't need to invoke actual parallelism, just show that the component-wise approach is correct
#[test]
fn watterson_theta_equivalent_concurrent(
    seed in 0..u64::MAX,
    ploidy in 1_usize..5,
    num_samples in 2_usize..50,
    missing_data_rate_raw in 0_f64..1.0 - f64::EPSILON,
    num_sites in 1_usize..10,
    non_normalized_freqs in vec(vec(f64::EPSILON..1_f64, 1..10), 10),
) {
    let mut rng = StdRng::seed_from_u64(seed);
    let sites = super::testdata::make_random_sites(
        &mut rng,
        ploidy,
        num_samples,
        missing_data_rate_raw,
        num_sites,
        non_normalized_freqs,
    );

    let counts = crate::testing::testdata::single_pop_counts(&mut sites.iter());
    let theta = counts.iter().try_fold(WattersonTheta::default(), |mut theta, s| {
        theta.try_add_site(s)?;
        Ok::<_, PopgenError>(theta)
    });

    let mut components = counts.iter().map(|s| {
        let mut theta = WattersonTheta::default();
        theta.try_add_site(s)?;
        Ok::<_, PopgenError>(theta)
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
        Ok(WattersonTheta::default())
    };

    match theta {
        Err(e) => {
            assert_eq!(std::mem::discriminant(&e), std::mem::discriminant(&parallel.err().unwrap()));
        },
        Ok(pi) => {
            let theta_raw = pi.as_raw();
            let parallel_raw = parallel.unwrap().as_raw();
            assert!((theta_raw - parallel_raw).abs() < 0.00001)
        },
    }
}
);
