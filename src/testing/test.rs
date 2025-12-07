mod tests {
    use crate::AlleleID;

    use crate::counts::MultiSiteCounts;
    use crate::iter::SiteCounts;
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

    // inefficient O(n^2) computation of pairwise diversity at a site
    //
    // hidden in test module because nobody should use this; just want to verify without magic numbers that calculations are correct
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
}
