use crate::stats::GlobalStatistic;
use crate::stats::TajimaD;
use crate::AlleleID;
use crate::MultiSiteCounts;
use rand::rng;
use rand::rngs::ThreadRng;
use rand::seq::SliceRandom;
use std::iter::repeat_n;

// NOTE: the copy/pase of shuffled_site can be
// considered okay b/c this test will go away
// once we have a naive implementation in place
// and can test using our random data API.
#[test]
fn tajima_d() {
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

    let allele_counts = MultiSiteCounts::from_tabular(sites);

    let tajima = TajimaD::from_iter_sites(allele_counts.iter());
    assert!((tajima.as_raw() - -0.15474069911037955).abs() < f64::EPSILON);
}
