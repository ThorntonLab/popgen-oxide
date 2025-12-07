use crate::AlleleID;

use crate::counts::MultiSiteCounts;
use crate::iter::SiteCounts;
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

    let counts = MultiSiteCounts::from_tabular(sites);

    assert_eq!(counts.len(), 2);
    assert!(!counts.is_empty());

    let mut iter = counts.iter();
    assert_eq!(
        iter.next().unwrap(),
        SiteCounts {
            counts: &[8, 7],
            total_alleles: 8 + 7 + 4,
        }
    );
    assert_eq!(
        iter.next().unwrap(),
        SiteCounts {
            counts: &[341, 69, 926],
            total_alleles: 341 + 69 + 926 + 300,
        }
    );
    assert!(iter.next().is_none());
}

#[test]
fn empty_counts() {
    let counts = MultiSiteCounts::default();

    assert!(counts.is_empty());
    assert_eq!(counts.len(), 0);
}

#[test]
#[should_panic]
fn bad_site_negative_count() {
    let mut counts = MultiSiteCounts::default();

    counts.add_site_from_counts([-1, -2, -3], 100).unwrap();
}

#[test]
#[should_panic]
fn bad_site_deficient_total() {
    let mut counts = MultiSiteCounts::default();

    counts.add_site_from_counts([1, 2, 3], 1).unwrap();
}
