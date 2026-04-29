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
fn pi_from_random_data(
    seed in 0..u64::MAX,
    ploidy in 1_usize..10,
    num_samples in 1_usize..50,
    missing_data_rate_raw in 0_f64..1.0 - f64::EPSILON,
    num_sites in 1_usize..10,
    non_normalized_freqs in vec(vec(f64::EPSILON..1_f64, 1..10), 10),
) {
    use rand::prelude::*;

    let mut rng = StdRng::seed_from_u64(seed);
    let mut sites = super::testdata::make_random_sites(
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
    let pi_from_counts = GlobalPi::try_from_iter_sites(counts.iter());
    let pi_naive = crate::testing::naivecalculations::pi(sites.iter());
    // compare
    match pi_from_counts {
        Err(_) => assert!(pi_naive.is_nan(), "{pi_naive} {counts:?}"),
        Ok(value) => assert!(
            (value.as_raw() - pi_naive).abs() <= 1e-10,
            "{pi_from_counts:?} != {pi_naive} {counts:?}"
        ),
    }
}
);

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
                let pi_from_counts = GlobalPi::try_from_iter_sites(counts.iter());
                assert_eq!(pi_from_counts.unwrap().as_raw(), 0.);
            }
        }
    }
}

#[test]
fn pi_try_from_iter_empty_is_err() {
    let c = crate::MultiSiteCounts::default();
    assert!(GlobalPi::try_from_iter_sites(c.iter()).is_err());
}
