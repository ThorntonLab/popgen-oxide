use crate::stats::GlobalStatistic;
use crate::stats::{Diversity, FStatistics};
use crate::testing::testdata::RandomSiteOptions;
use crate::{MultiPopulationCounts, PopgenError};
use std::borrow::Cow;
use std::collections::HashMap;

#[test]
fn f_st_no_pops() {
    // equivalently, could have populations but exclude them all...
    // but this is simpler
    let counts = MultiPopulationCounts::default();
    assert!(matches!(
        FStatistics::try_from_populations(&counts, |_| Some(1.0)),
        Err(PopgenError::CalculationError)
    ));
}

#[test]
fn f_st_empty_pops() {
    for n_pops in [1, 2, 5] {
        assert!(matches!(
            FStatistics::try_from_populations(
                &MultiPopulationCounts::of_empty_populations(n_pops),
                |_| { Some(1.0) }
            ),
            Err(PopgenError::EmptySiteCounts)
        ));
    }
}

#[test]
fn f_st() {
    let mut populations = MultiPopulationCounts::of_empty_populations(3);

    let data = [([1, 2, 0], 3), ([3, 0, 0], 3), ([0, 1, 2], 3)];
    let weights = [1.0, 2.0, 3.0];

    populations
        .extend_populations_from_site(|i| (&data[i].0, data[i].1))
        .unwrap();

    let f_st = FStatistics::try_from_populations(&populations, |i| Some(weights[i])).unwrap();

    let (pi_b_top, pi_b_bottom) = f_st.pi_b_parts();

    assert!(
        (pi_b_top
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
        pi_b_bottom,
        weights
            .iter()
            .enumerate()
            .flat_map(|(i, w1)| weights.iter().skip(i + 1).map(move |w2| w1 * w2))
            .sum::<f64>()
    );

    let (pi_s_top, pi_s_bottom) = f_st.pi_s_parts();

    assert!(
        (pi_s_top
            - (0..populations.num_populations())
                .map(|pop_i| {
                    // sum of weight * weight * diversity within this population
                    Diversity::try_from_iter_sites(populations.iter_sites_in(pop_i))
                        .unwrap()
                        .as_raw()
                        * weights[pop_i].powi(2)
                })
                .sum::<f64>())
        .abs()
            < f64::EPSILON
    );

    assert!(
        (pi_s_bottom
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
                        (Cow::Owned(counts_vec), alleles.len() as i32)
                    })
                    .unwrap();
            }

            let f_st_from_counts =
                FStatistics::try_from_populations(&counts, |i| Some(pop_weights[i])).unwrap();
            let (pi_total_naive, pi_self_naive, pi_between_naive) =
                crate::testing::naivecalculations::f_st(
                    &mut pops
                        .iter()
                        .cloned()
                        .zip(pop_weights.iter())
                        .map(|(p, &w)| (w, p)),
                );
            assert!((1.0 - (f_st_from_counts.pi_t() / pi_total_naive)).abs() < 0.00001);
            assert!((1.0 - (f_st_from_counts.pi_s().unwrap() / pi_self_naive)).abs() < 0.00001);
            assert!((1.0 - (f_st_from_counts.pi_b().unwrap() / pi_between_naive)).abs() < 0.00001);
            for i in 0..n_pops {
                for j in i + 1..n_pops {
                    let naivef2 =
                        crate::testing::naivecalculations::f2(i, j, &mut pops.iter().cloned());
                    let f2 = f_st_from_counts.f2(i, j).unwrap();
                    assert!((naivef2 - f2).abs() < 1e-6, "{i} {j}: {naivef2} != {f2}");
                }
            }
        }
    }
}
