use crate::stats::StatRepresentation;
use crate::stats::TajimasD;
use crate::stats::UnpolarisedSiteStat;
use crate::AlleleID;
use crate::SampleAlleleCounts;
use rand::rng;
use rand::rngs::ThreadRng;
use rand::seq::SliceRandom;
use std::iter::repeat_n;

use proptest::collection::vec;
use proptest::prelude::*;

// NOTE: the copy/pase of shuffled_site can be
// considered okay b/c this test will go away
// once we have a naive implementation in place
// and can test using our random data API.
#[test]
fn tajimas_d() {
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
        // common mutation
        vec![(Some(AlleleID::from(0)), 11), (Some(AlleleID::from(1)), 7)],
        // rare mutation
        vec![(Some(AlleleID::from(0)), 16), (Some(AlleleID::from(1)), 2)],
        // rare mutation
        vec![(Some(AlleleID::from(0)), 1), (Some(AlleleID::from(1)), 17)],
    ]
    .into_iter()
    .map(|site| shuffled_site(site.into_iter(), &mut rng))
    .collect::<Vec<_>>();

    let allele_counts = SampleAlleleCounts::try_from_tabular(sites).unwrap();

    let d = TajimasD::try_from_iter_sites(allele_counts.iter()).unwrap();
    assert!((d.as_raw() - -0.15474069911037955).abs() < f64::EPSILON);
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
        let diversity_from_counts = TajimasD::try_from_iter_sites(counts.iter());
        if let Ok(value) = diversity_from_counts {
            if !value.as_raw().is_nan() {
                let splitlen = counts.len() / max_num_splits;
                let mut div_split = vec![];
                for nsplits in 0..max_num_splits {
                    let takelen = if nsplits < max_num_splits - 1 {
                        splitlen
                    } else {
                        counts.len() - nsplits * splitlen
                    };
                    let div = TajimasD::try_from_iter_sites(counts.iter().skip(nsplits * splitlen).take(takelen)).unwrap_or_default();
                    div_split.push(div);
                }
                let reduced = div_split.iter().fold(TajimasD::default(), |acc, &i| acc.try_reduce(i).unwrap());
                let absdiff = (value.as_raw() - reduced.as_raw()).abs();
                assert!(absdiff <= 1e-9, "{value:?} != {reduced:?} ({} {} {absdiff}) ({div_split:?}) {splitlen} {counts:?}", value.as_raw(), reduced.as_raw())
            }
        }
    }
);
