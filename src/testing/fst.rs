use crate::stats::FStatisticParts;
use crate::stats::GlobalPi;
use crate::stats::GlobalStatistic;
use crate::MultiPopulationCounts;

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
                    GlobalPi::from_iter_sites(pop.iter()).as_raw() * (weights[i]).powi(2)
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
