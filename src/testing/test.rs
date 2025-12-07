mod tests {
    use crate::{AlleleID, Count};

    use crate::counts::{MultiPopulationCounts, MultiSiteCounts};
    use crate::iter::SiteCounts;
    use crate::stats::{FStatisticParts, GlobalPi, GlobalStatistic, TajimaD, WattersonTheta};
    use itertools::Itertools;
    use rand::rng;
    use rand::rngs::ThreadRng;
    use rand::seq::SliceRandom;
    use std::iter::repeat_n;
    // use triangle_matrix::{
    //     SimpleLowerTri, SymmetricUpperTri, SymmetricUpperTriMut, Triangle, TriangleMut,
    // };

    // struct TriVec<T>(usize, Vec<T>);

    // impl<T> Triangle<T> for TriVec<T> {
    //     type Inner = Vec<T>;

    //     fn n(&self) -> usize {
    //         self.0
    //     }

    //     fn inner(&self) -> &Self::Inner {
    //         &self.1
    //     }
    // }

    // impl<T> TriangleMut<T> for TriVec<T> {
    //     fn inner_mut(&mut self) -> &mut Self::Inner {
    //         &mut self.1
    //     }
    // }

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

    /// inefficient O(n^2) computation of pairwise diversity at a site
    ///
    /// hidden in test module because nobody should use this; just want to verify without magic numbers that calculations are correct
    // fn pi_from_matrix(alleles: &[Option<AlleleID>]) -> f64 {
    //     use SymmetricUpperTri;
    //     use SymmetricUpperTriMut;

    //     let total_alleles = alleles.len();
    //     let mut mat = TriVec(
    //         total_alleles,
    //         Vec::from_iter(repeat_n(
    //             None::<i32>,
    //             total_alleles * (total_alleles + 1) / 2,
    //         )),
    //     );

    //     for (ind_a, allele_a) in alleles.iter().enumerate() {
    //         for (ind_b, allele_b) in alleles.iter().enumerate().skip(ind_a + 1) {
    //             *SymmetricUpperTriMut::get_element_mut(&mut mat, ind_a, ind_b) = match *allele_a {
    //                 None => None,
    //                 Some(allele_id_a) => match *allele_b {
    //                     None => None,
    //                     Some(allele_id_b) => match allele_id_a.eq(&allele_id_b) {
    //                         false => Some(1),
    //                         true => Some(0),
    //                     },
    //                 },
    //             }
    //         }
    //     }

    //     let make_iter = || {
    //         mat.iter_triangle_indices()
    //             .filter_map(|(i, j)| *SymmetricUpperTri::get_element(&mat, i, j))
    //     };
    //     make_iter().sum::<i32>() as f64 / make_iter().count() as f64
    // }

    #[test]
    fn load_raw() {
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

    // #[test]
    // fn global_pi() {
    //     let mut rng = rng();
    //     let sites = vec![
    //         shuffled_site(
    //             vec![
    //                 (Some(AlleleID::from(0)), 35),
    //                 (Some(AlleleID::from(1)), 6),
    //                 (None, 3),
    //             ]
    //             .into_iter(),
    //             &mut rng,
    //         ),
    //         shuffled_site(
    //             vec![
    //                 (Some(AlleleID::from(0)), 2),
    //                 (Some(AlleleID::from(1)), 14),
    //                 (Some(AlleleID::from(2)), 155),
    //             ]
    //             .into_iter(),
    //             &mut rng,
    //         ),
    //     ];

    //     let expect_site_0 = pi_from_matrix(&sites[0]);
    //     let expect_site_1 = pi_from_matrix(&sites[1]);

    //     let allele_counts = MultiSiteCounts::from_tabular(sites);

    //     assert!(
    //         (GlobalPi::from_iter_sites(allele_counts.iter().take(1)).as_raw() - expect_site_0)
    //             .abs()
    //             < f64::EPSILON
    //     );
    //     assert!(
    //         (GlobalPi::from_iter_sites(allele_counts.iter().skip(1).take(1)).as_raw()
    //             - expect_site_1)
    //             .abs()
    //             < f64::EPSILON
    //     );
    //     assert!(
    //         (GlobalPi::from_iter_sites(allele_counts.iter()).as_raw()
    //             - (expect_site_0 + expect_site_1))
    //             .abs()
    //             < f64::EPSILON
    //     );
    // }

    #[test]
    fn tajima_d() {
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
        .collect_vec();

        let allele_counts = MultiSiteCounts::from_tabular(sites);

        let tajima = TajimaD::from_iter_sites(allele_counts.iter());
        assert!((tajima.as_raw() - -0.15474069911037955).abs() < f64::EPSILON);
    }

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
}
