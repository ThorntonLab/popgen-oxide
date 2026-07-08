use crate::stats::StatRepresentation;
use crate::stats::UnpolarisedSiteStat;
use crate::stats::WattersonsTheta;

use proptest::collection::vec;
use proptest::prelude::*;

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
        let theta = WattersonsTheta::try_from_iter_sites(counts.iter());
        let theta_naive = crate::testing::naivecalculations::watterson_theta(&mut sites.iter_mut());
        match theta {
            Err(_) => assert!(theta_naive.is_nan()),
            Ok(value) => assert!(
                (value.as_raw() - theta_naive).abs() <= 1e-10,
                "{value:?} != {theta_naive}"
            ),
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
            let theta = WattersonsTheta::try_from_iter_sites(counts.iter());
            let theta_naive =
                crate::testing::naivecalculations::watterson_theta(&mut sites.iter_mut());
            // compare
            match theta {
                Err(_) => assert!(theta_naive.is_nan(), "{theta_naive}"),
                Ok(value) => assert!(
                    (value.as_raw() - theta_naive).abs() <= 1e-10,
                    "{value:?} != {theta_naive}"
                ),
            }
        }
    }
}

#[test]
fn wattherson_theta_try_from_iter_empty_is_err() {
    let c = crate::SampleAlleleCounts::default();
    assert!(WattersonsTheta::try_from_iter_sites(c.iter()).is_err());
}

proptest!(
    #[test]
    fn test_try_reduce(
        seed in 0..u64::MAX,
        ploidy in 1_usize..10,
        num_samples in 1_usize..50,
        missing_data_rate_raw in 0_f64..1.0 - f64::EPSILON,
        num_sites in 4_usize..10,
        non_normalized_freqs in vec(vec(f64::EPSILON..1_f64, 1..10), 10),
        max_num_splits in 1_usize..4
    ) {
        use rand::prelude::*;
        use crate::traits::TryReduce;

        let mut rng = StdRng::seed_from_u64(seed);
        let sites = super::testdata::make_random_sites(
            &mut rng,
            ploidy,
            num_samples,
            missing_data_rate_raw,
            num_sites,
            non_normalized_freqs,
        );
        // convert to our normal format
        let counts = crate::testing::testdata::single_pop_counts(&mut sites.iter());

        // get the calcs
        let diversity_from_counts = WattersonsTheta::try_from_iter_sites(counts.iter());
        if let Ok(value) = diversity_from_counts {
            let splitlen = counts.len() / max_num_splits;
            let mut div_split = vec![];
            for nsplits in 0..max_num_splits {
                let takelen = if nsplits < max_num_splits - 1 {
                    splitlen
                } else {
                    counts.len() - nsplits * splitlen
                };
                let div = WattersonsTheta::try_from_iter_sites(counts.iter().skip(nsplits * splitlen).take(takelen)).unwrap_or_default();
                div_split.push(div);
            }
            let reduced = div_split.iter().fold(WattersonsTheta::default(), |acc, &i| acc.try_reduce(i).unwrap());
            assert!((value.as_raw() - reduced.as_raw()).abs() <= 1e-9, "{value:?} != {reduced:?} ({div_split:?}) {splitlen} {counts:?}")
        }
    }
);
