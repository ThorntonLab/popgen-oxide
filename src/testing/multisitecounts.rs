use crate::{AlleleCounts, AlleleID};

use crate::counts::SampleAlleleCounts;
use proptest::collection::vec;
use proptest::prelude::*;
use rand::rng;
use rand::rngs::ThreadRng;
use rand::seq::SliceRandom;
use std::iter::repeat_n;

#[test]
fn load_raw() {
    // NOTE: copy-paste of a fn present elsewhere
    fn shuffled_site(
        ids: impl Iterator<Item = (Option<AlleleID>, usize)>,
        rng: &mut ThreadRng,
    ) -> Vec<Option<AlleleID>> {
        let mut site = vec![];
        ids.for_each(|(id, count)| {
            site.extend(repeat_n(id, count));
        });

        site.shuffle(rng);
        site
    }
    let mut rng = rng();

    let sites = vec![
        shuffled_site(
            vec![(Some(AlleleID(0)), 8), (Some(AlleleID(1)), 7), (None, 4)].into_iter(),
            &mut rng,
        ),
        shuffled_site(
            vec![
                (Some(AlleleID(0)), 341),
                (Some(AlleleID::from(1)), 69),
                (Some(AlleleID::from(2)), 926),
                (None, 300),
            ]
            .into_iter(),
            &mut rng,
        ),
    ];

    let counts = SampleAlleleCounts::try_from_tabular(sites).unwrap();

    assert_eq!(counts.len(), 2);
    assert!(!counts.is_empty());

    let mut iter = counts.iter();
    assert_eq!(
        iter.next().unwrap(),
        AlleleCounts {
            counts: &[8, 7],
            total_alleles: 8 + 7 + 4,
        }
    );
    assert_eq!(
        iter.next().unwrap(),
        AlleleCounts {
            counts: &[341, 69, 926],
            total_alleles: 341 + 69 + 926 + 300,
        }
    );
    assert!(iter.next().is_none());
}

#[test]
fn empty_counts() {
    let counts = SampleAlleleCounts::default();

    assert!(counts.is_empty());
    assert_eq!(counts.len(), 0);
}

#[test]
#[should_panic]
fn bad_site_negative_count() {
    AlleleCounts::try_new(&[-1, -2, -3], 100).unwrap();
}

#[test]
#[should_panic]
fn bad_site_empty_count() {
    AlleleCounts::try_new(&[], 100).unwrap();
}

#[test]
#[should_panic]
fn bad_site_deficient_total() {
    AlleleCounts::try_new(&[1, 2, 3], 1).unwrap();
}

fn test_try_reduce_details(
    seed: u64,
    ploidy: usize,
    num_samples: usize,
    missing_data_rate_raw: f64,
    num_sites: usize,
    non_normalized_freqs: Vec<Vec<f64>>,
    max_num_splits: usize,
) {
    use crate::traits::TryReduce;
    use rand::prelude::*;

    let mut rng = StdRng::seed_from_u64(seed);
    let mut splitcounts = vec![];
    for _ in 0..max_num_splits + 1 {
        let sites = super::testdata::make_random_sites(
            &mut rng,
            ploidy,
            num_samples,
            missing_data_rate_raw,
            num_sites,
            non_normalized_freqs.clone(),
        );
        // convert to our normal format
        let counts = crate::testing::testdata::single_pop_counts(&mut sites.iter());
        splitcounts.push(counts);
    }

    let mut mergedcounts = crate::counts::SampleAlleleCounts::default();
    for i in splitcounts.iter().flat_map(|c| c.iter()) {
        mergedcounts.add_site_from_counts(i);
    }
    let reduced_counts = splitcounts
        .into_iter()
        .try_fold(crate::counts::SampleAlleleCounts::default(), |a, b| {
            a.try_reduce(b)
        })
        .unwrap();
    assert_eq!(reduced_counts.len(), mergedcounts.len());
    for (i, j) in reduced_counts.iter().zip(mergedcounts.iter()) {
        assert_eq!(i.counts(), j.counts());
        assert_eq!(i.total_alleles(), j.total_alleles());
    }
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
        test_try_reduce_details(seed, ploidy, num_samples, missing_data_rate_raw, num_sites, non_normalized_freqs, max_num_splits);
    }
);
