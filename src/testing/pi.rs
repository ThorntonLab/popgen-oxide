use crate::stats::GlobalPi;
use crate::stats::GlobalStatistic;

use proptest::collection::vec;
use proptest::prelude::*;

proptest!(
// The basic logic here is that we should
// be able to make random data, convert it
// to our normal format, and get stats
// from both data inputs and get the same
// answer.
// If this is not possible then we have bugs
// in one of several possible places...
#[test]
fn pi_from_random_data(seed in 0..u64::MAX,
                       ploidy in 1_usize..10,
                       num_samples in 1_usize..50,
                       missing_data_rate_raw in 0_f64..1.0 - f64::EPSILON,
                       non_normalized_freqs in vec(0_f64..1_f64, 1..10)) {
    use rand::prelude::*;

    let mut rng = StdRng::seed_from_u64(seed);
    let sum = non_normalized_freqs.iter().sum::<f64>();
    let freqs = non_normalized_freqs.into_iter().map(|fi|fi/sum).collect::<Vec<_>>();
    let mut sites = vec![];
    // make num_samples random sites.
    // No missing data, etc..
    let missing_data_rate = if missing_data_rate_raw > 0.0 {
        Some(crate::testing::testdata::RandomSiteOptions{missing_data_rate:Some(missing_data_rate_raw)})
    } else { None};
    let site =
        crate::testing::testdata::random_site_rng(num_samples, ploidy, &freqs, missing_data_rate, &mut rng);
    sites.push(site);
    // convert to our normal format
    let counts = crate::testing::testdata::single_pop_counts(&mut sites.iter());
    // get the calcs
    let pi_from_counts = GlobalPi::from_iter_sites(counts.iter());
    let pi_naive = crate::testing::naivecalculations::pi(&mut sites.iter_mut());
    // compare
    if pi_naive.is_nan() {
        assert!(pi_from_counts.as_raw().is_nan(), "{pi_from_counts:?} {counts:?}");
    } else {
        assert!(
            (pi_from_counts.as_raw() - pi_naive).abs() <= 1e-10,
            "{pi_from_counts:?} != {pi_naive}"
        );
    }
}
);

#[test]
fn pi_from_random_data_with_missing_data() {
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
            let pi_from_counts = GlobalPi::from_iter_sites(counts.iter());
            let pi_naive = crate::testing::naivecalculations::pi(&mut sites.iter_mut());
            // compare
            if pi_naive.is_nan() {
                assert!(pi_from_counts.as_raw().is_nan());
            } else {
                assert!(
                    (pi_from_counts.as_raw() - pi_naive).abs() <= 1e-10,
                    "{pi_from_counts:?} != {pi_naive}"
                );
            }
        }
    }
}

#[test]
fn pi_allele_frequency_of_one() {
    use rand::prelude::*;

    let mut rng = StdRng::seed_from_u64(54321);
    for ploidy in [1, 2, 4] {
        for num_alleles in [2, 3, 4] {
            let freqs_iter = crate::testing::testdata::FixedMutationIterator::new(num_alleles);

            for freqs in freqs_iter.iter() {
                assert_eq!(freqs.len(), num_alleles);
                let site =
                    crate::testing::testdata::random_site_rng(10, ploidy, &freqs, None, &mut rng);
                // convert to our normal format
                let counts =
                    crate::testing::testdata::single_pop_counts(&mut std::iter::once(&site));
                let pi_from_counts = GlobalPi::from_iter_sites(counts.iter());
                assert_eq!(pi_from_counts.as_raw(), 0.);
            }
        }
    }
}
