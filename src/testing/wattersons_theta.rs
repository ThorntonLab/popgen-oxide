use crate::stats::GlobalStatistic;
use crate::stats::WattersonTheta;
use crate::AlleleID;
use crate::Count;
use crate::MultiSiteCounts;

#[test]
fn watterson_theta() {
    let sites = vec![
        vec![Some(0), Some(0), Some(1), Some(1), Some(1), Some(2)],
        vec![Some(0), Some(1), Some(1), Some(1), Some(2), None],
    ]
    .into_iter()
    .map(|site| {
        site.into_iter()
            .map(|sam| sam.map(AlleleID::from))
            .collect::<Vec<_>>()
    })
    .collect::<Vec<_>>();

    let allele_counts = MultiSiteCounts::from_tabular(sites);
    let theta = WattersonTheta::from_iter_sites(allele_counts.iter());

    let mut expected = 0f64;
    for site in allele_counts.iter() {
        let num_variants = site.counts.iter().filter(|c| **c > 0).count();
        let total_samples = site.counts.iter().filter(|c| **c > 0).sum::<Count>();

        if num_variants > 1 {
            let numerator = (num_variants - 1) as f64;
            let denominator = (1..total_samples).map(|i| 1f64 / i as f64).sum::<f64>();
            expected += numerator / denominator;
        }
    }

    assert!((theta.as_raw() - expected).abs() < f64::EPSILON);
}

#[test]
fn watterson_theta_from_random_data() {
    use rand::prelude::*;

    let mut rng = StdRng::seed_from_u64(54321);
    for ploidy in [1, 2, 3, 4] {
        let freqs = [0.25, 0.5, 0.25]; // fixed allele freqs per site
        let mut sites = vec![];
        // make 10 random sites.
        // No missing data, etc..
        for _ in 0..10 {
            let site =
                crate::testing::testdata::random_site_rng(10, ploidy, &freqs, None, &mut rng);
            sites.push(site);
        }
        // convert to our normal format
        let counts = crate::testing::testdata::single_pop_counts(&mut sites.iter());
        let theta = WattersonTheta::from_iter_sites(counts.iter());
        let theta_naive = crate::testing::naivecalculations::watterson_theta(&mut sites.iter_mut());
        if theta_naive.is_nan() {
            assert!(theta.as_raw().is_nan())
        } else {
            assert!((theta.as_raw() - theta_naive).abs() <= 1e-10)
        }
    }
}

#[test]
fn watterson_theta_from_random_data_with_missing_data() {
    use rand::prelude::*;

    let mut rng = StdRng::seed_from_u64(54321);

    for ploidy in [1, 2, 4] {
        for rate in [0.01, 0.1, 0.5, 0.9] {
            let freqs = [0.25, 0.5, 0.25]; // fixed allele freqs per site
            let mut sites = vec![];
            // make 10 random sites.
            // No missing data, etc..
            for _ in 0..10 {
                let site = crate::testing::testdata::random_site_rng(
                    10,
                    ploidy,
                    &freqs,
                    Some(crate::testing::testdata::RandomSiteOptions {
                        missing_data_rate: Some(rate),
                    }),
                    &mut rng,
                );
                sites.push(site);
            }
            // convert to our normal format
            let counts = crate::testing::testdata::single_pop_counts(&mut sites.iter());
            // get the calcs
            let theta = WattersonTheta::from_iter_sites(counts.iter());
            let theta_naive =
                crate::testing::naivecalculations::watterson_theta(&mut sites.iter_mut());
            // compare
            if theta_naive.is_nan() {
                assert!(theta.as_raw().is_nan());
            } else {
                assert!(
                    (theta.as_raw() - theta_naive).abs() <= 1e-10,
                    "{theta:?} != {theta_naive}"
                );
            }
        }
    }
}
