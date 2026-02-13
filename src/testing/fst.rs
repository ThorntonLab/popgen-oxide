use crate::stats::FStatisticParts;
use crate::stats::GlobalPi;
use crate::stats::GlobalStatistic;
use crate::testing::testdata::RandomSiteOptions;
use crate::MultiPopulationCounts;
use std::borrow::Cow;
use std::collections::HashMap;

#[test]
fn f_st_empty() {
    fn ok(populations: &MultiPopulationCounts) {
        // this is the one case where these fail; let's make sure that is the case
        let f_st = populations.f_st_if(|_| None);
        assert_eq!(f_st.pi_S(), None);
        assert_eq!(f_st.pi_B(), None);
        assert_eq!(f_st.pi_D(), None);
    }

    ok(&MultiPopulationCounts::of_empty_populations(0));
    ok(&MultiPopulationCounts::of_empty_populations(1));
    ok(&MultiPopulationCounts::of_empty_populations(5));
}

#[test]
fn f_st() {
    let mut populations = MultiPopulationCounts::of_empty_populations(3);

    let data = [([1, 2, 0], 3), ([3, 0, 0], 3), ([0, 1, 2], 3)];
    let weights = [1.0, 2.0, 3.0];

    populations
        .extend_populations_from_site(|i| (&data[i].0, data[i].1))
        .unwrap();

    let f_st = populations.f_st_if(|i| Some(weights[i]));

    #[allow(non_snake_case)]
    let (pi_B_top, pi_B_bottom) = f_st.pi_B_parts();

    assert!(
        (pi_B_top
                - (
                // differences between (0, 1)
                6.0 * weights[0] * weights[1]
                    // (0, 2)
                    + 7.0 * weights[0] * weights[2]
                    // (1, 2)
                    + 9.0 * weights[1] * weights[2]
            )
                // number of comparisons between any two populations
                / 9.)
                .abs()
                // this comparison seems a little finicky
                < 10.0_f64.powi(-10)
    );

    assert_eq!(
        pi_B_bottom,
        weights
            .iter()
            .enumerate()
            .flat_map(|(i, w1)| weights.iter().skip(i + 1).map(move |w2| w1 * w2))
            .sum::<f64>()
    );

    #[allow(non_snake_case)]
    let (pi_S_top, pi_S_bottom) = f_st.pi_S_parts();

    assert!(
        (pi_S_top
            - populations
                .iter()
                .enumerate()
                .map(|(i, pop)| {
                    // sum of weight * weight * pi within this population
                    GlobalPi::from_iter_sites(pop.iter()).unwrap().as_raw() * (weights[i]).powi(2)
                })
                .sum::<f64>())
        .abs()
            < f64::EPSILON
    );

    assert!(
        (pi_S_bottom
                // denominator should be sum of squares of weights
                - weights.iter().map(|w| w.powi(2)).sum::<f64>())
        .abs()
            < f64::EPSILON
    );
}

#[test]
fn f_st_from_random_data() {
    use rand::prelude::*;
    let n_sites = 5;
    let pop_weights = [0.1, 0.2, 0.3, 0.4];
    let n_alleles = 3;

    let mut rng = StdRng::seed_from_u64(54321);
    for n_pops in 2..4 {
        for ploidy in 1..3 {
            let pops = (0..n_pops)
                .map(|_| {
                    let mut sites = vec![];
                    for _ in 0..n_sites {
                        let site = crate::testing::testdata::random_site_rng(
                            5,
                            ploidy,
                            &[0.25, 0.5, 0.25],
                            Some(RandomSiteOptions {
                                missing_data_rate: Some(0.1),
                            }),
                            &mut rng,
                        );
                        sites.push(site);
                    }
                    sites
                })
                .collect::<Vec<_>>();

            let mut counts = MultiPopulationCounts::of_empty_populations(n_pops);
            #[expect(
                clippy::needless_range_loop,
                reason = "https://github.com/rust-lang/rust-clippy/issues/16344"
            )]
            for s in 0..n_sites {
                counts
                    .extend_populations_from_site(|pop_i| {
                        let alleles = pops[pop_i][s]
                            .iter()
                            .flat_map(|gts| gts.iter())
                            .collect::<Vec<_>>();
                        let mut counts = HashMap::new();
                        for non_missing_allele_id in alleles.iter().flatten() {
                            *counts.entry(non_missing_allele_id).or_default() += 1;
                        }

                        let counts_vec = (0..n_alleles)
                            .map(|allele_id| counts.get(&allele_id).cloned().unwrap_or_default())
                            .collect::<Vec<_>>();
                        (Cow::Owned(counts_vec), alleles.len())
                    })
                    .unwrap();
            }

            let f_st_from_counts = counts.f_st_if(|i| Some(pop_weights[i]));
            let (pi_total_naive, pi_self_naive, pi_between_naive) =
                crate::testing::naivecalculations::f_st(
                    &mut pops
                        .iter()
                        .cloned()
                        .zip(pop_weights.iter())
                        .map(|(p, &w)| (w, p)),
                );
            assert!((1.0 - (f_st_from_counts.pi_T() / pi_total_naive)).abs() < 0.00001);
            assert!((1.0 - (f_st_from_counts.pi_S().unwrap() / pi_self_naive)).abs() < 0.00001);
            assert!((1.0 - (f_st_from_counts.pi_B().unwrap() / pi_between_naive)).abs() < 0.00001);
        }
    }
}
