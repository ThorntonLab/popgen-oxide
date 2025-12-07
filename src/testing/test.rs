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

    #[cfg(feature = "noodles")]
    mod vcf {
        use crate::adapter::vcf::{
            record_to_genotypes_adapter, VCFToPopulationsAdapter, WhichPopulation,
        };
        use crate::counts::MultiSiteCounts;
        use crate::{AlleleID, PopgenResult};
        use noodles::vcf::header::record::value::map::{Contig, Format};
        use noodles::vcf::header::record::value::Map;
        use noodles::vcf::variant::io::Write;
        use noodles::vcf::variant::record::samples::keys::key;
        use noodles::vcf::variant::record::samples::series::value::genotype::Phasing::Unphased;
        use noodles::vcf::variant::record_buf::samples::sample::value::genotype::Allele;
        use noodles::vcf::variant::record_buf::samples::sample::value::Genotype;
        use noodles::vcf::variant::record_buf::samples::sample::Value;
        use noodles::vcf::variant::record_buf::samples::Keys;
        use noodles::vcf::variant::record_buf::{AlternateBases, Samples};
        use noodles::vcf::variant::RecordBuf;
        use rand::prelude::SliceRandom;
        use rand::rng;
        use std::borrow::Cow;
        use std::collections::HashMap;

        // SiteVariant is to be a slice like ["A", "AG"] for a sample with these two genotypes
        // the appropriate IDs, number of samples, etc. will be calculated
        #[allow(dead_code)]
        fn make_record<SiteGenotypes, SiteVariant, InnerAllele>(
            seq_name: &str,
            site_genotypes: SiteGenotypes,
        ) -> Option<RecordBuf>
        where
            SiteGenotypes: AsRef<[(SiteVariant, usize)]>,
            SiteVariant: AsRef<[Option<InnerAllele>]>,
            InnerAllele: AsRef<str>,
        {
            let site_genotypes = site_genotypes.as_ref();
            if site_genotypes.len() < 2 {
                return None;
            }

            let mut alleles_seen = Vec::new();
            for (genotype, _) in site_genotypes {
                for gt_str in genotype
                    .as_ref()
                    .iter()
                    .filter_map(|gt_str| gt_str.as_ref())
                    .map(|gt_str| gt_str.as_ref())
                {
                    if !alleles_seen.iter().any(|allele| gt_str.eq(*allele)) {
                        alleles_seen.push(gt_str);
                    }
                }
            }

            Some(
                RecordBuf::builder()
                    .set_reference_sequence_name(seq_name)
                    .set_variant_start(noodles_core::Position::MIN)
                    .set_reference_bases(alleles_seen[0])
                    .set_alternate_bases(AlternateBases::from(
                        alleles_seen[1..]
                            .iter()
                            .map(|s| String::from(*s))
                            .collect::<Vec<_>>(),
                    ))
                    .set_samples(Samples::new(
                        Keys::from_iter(vec![String::from(key::GENOTYPE)]),
                        // for each genotype and its count, create that many samples accordingly
                        {
                            let mut all_samples = site_genotypes
                                .iter()
                                .flat_map(|(genotype, count)| {
                                    vec![
                                        // for each sample at this site, there will only be the GT field and no other hence another layer of Vec here
                                        vec![Some(Value::from(Genotype::from_iter(
                                            genotype.as_ref().iter()
                                                .map(|sample_variant| Allele::new(
                                                    sample_variant.as_ref().map(|some| alleles_seen.iter()
                                                        .enumerate()
                                                        .find(|(_, variant)| **variant == some.as_ref())
                                                        .unwrap().0
                                                    ),
                                                    Unphased  // TODO: make this configurable?
                                                )),
                                        )))];
                                        *count
                                    ]
                                })
                                .collect::<Vec<_>>();
                            all_samples.shuffle(&mut rng());
                            all_samples
                        },
                    ))
                    .build(),
            )
        }

        #[allow(dead_code)]
        fn make_mock_vcf<Sites, SiteGenotypes, SiteVariant, Allele>(sites: Sites) -> Option<String>
        where
            Sites: AsRef<[SiteGenotypes]>,
            SiteGenotypes: AsRef<[(SiteVariant, usize)]>,
            SiteVariant: AsRef<[Option<Allele>]>,
            Allele: AsRef<str>,
        {
            let mut buf = Vec::new();
            let mut writer =
                noodles::vcf::io::writer::Builder::default().build_from_writer(&mut buf);

            let seq_name = "chr0";
            // how many samples are needed to contain the genotypes of the sample with the most genotypes?
            let num_samples = sites
                .as_ref()
                .iter()
                .map(|site_alleles| site_alleles.as_ref().iter().map(|(_, count)| *count).sum())
                .max()?;

            let header = {
                let id = key::GENOTYPE;
                let format = Map::<Format>::from(id);

                let mut header = noodles::vcf::Header::builder()
                    .add_contig(seq_name, Map::<Contig>::new())
                    .add_format(id, format.clone());

                for num in 0..num_samples {
                    header = header.add_sample_name(format!("s{num}"));
                }

                header.build()
            };

            writer.write_header(&header).unwrap();

            for site in sites.as_ref() {
                let record = make_record(seq_name, site).unwrap();
                writer.write_variant_record(&header, &record).unwrap();
            }

            // kick the writer out of the buffer so we can move it
            drop(writer);

            Some(String::from_utf8(buf).unwrap())
        }

        fn counts_from_vcf(
            vcf_buf: &str,
            ploidy: usize,
        ) -> (Vec<Vec<Option<AlleleID>>>, MultiSiteCounts) {
            let mut reader = noodles::vcf::io::reader::Builder::default()
                .build_from_reader(vcf_buf.as_bytes())
                .unwrap();

            let header = reader.read_header().unwrap();

            let all_alleles = reader
                .records()
                .map(Result::unwrap)
                .map(|rec| record_to_genotypes_adapter(&header, rec, ploidy))
                .collect::<PopgenResult<Vec<_>>>()
                .unwrap();
            let counts = MultiSiteCounts::from_tabular(all_alleles.iter().cloned());
            (all_alleles, counts)
        }

        fn make_vcf() -> &'static str {
            // let vcf_buf = make_mock_vcf(vec![
            //     vec![
            //         (vec![Some("A")], 11),
            //         (vec![Some("C")], 7),
            //     ],
            //     vec![
            //         (vec![Some("G")], 7),
            //         (vec![None], 3),
            //         (vec![Some("A")], 8),
            //     ]
            // ]).unwrap();
            //
            // println!("{}", &vcf_buf);

            r#"##fileformat=VCFv4.5
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##contig=<ID=chr0>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	s0	s1	s2	s3	s4	s5	s6	s7	s8	s9	s10	s11	s12	s13	s14	s15	s16	s17
chr0	1	.	A	C	.	.	.	GT	/0	/1	/1	/0	/1	/0	/1	/0	/0	/0	/0	/0	/0	/1	/0	/1	/1	/0
chr0	1	.	G	A	.	.	.	GT	/0	/1	/1	/0	/1	/1	/0	/0	/.	/.	/0	/0	/1	/1	/1	/1	/0	/."#
        }

        #[test]
        fn load_vcf() {
            let (_all_alleles, allele_counts) = counts_from_vcf(make_vcf(), 1);
            let mut iter = allele_counts.iter();
            let counts_0 = iter.next().unwrap();
            let counts_1 = iter.next().unwrap();
            assert!(iter.next().is_none());

            assert_eq!(counts_0.counts(), &[11, 7]);
            assert_eq!(counts_0.total_alleles(), 18);
            assert_eq!(counts_1.counts(), &[7, 8]);
            assert_eq!(counts_1.total_alleles(), 18);
        }

        #[test]
        fn load_vcf_multi_population() {
            // use this overly verbose and inefficient map to verify lifetime correctness
            let map = (0..18)
                .map(|n| (format!("s{}", n), ["A", "B"][n % 2].to_owned()))
                .collect::<HashMap<_, _>>();

            struct MapBasedMapper(HashMap<String, String>);

            impl WhichPopulation<()> for MapBasedMapper {
                fn which_population<'s>(
                    &'s mut self,
                    sample_name: &'_ str,
                ) -> Result<Cow<'s, str>, ()> {
                    self.0
                        .get(sample_name)
                        .map(String::as_str)
                        .map(Cow::Borrowed)
                        .ok_or(())
                }
            }

            let mut vcf_reader = noodles::vcf::io::reader::Builder::default()
                .build_from_reader(make_vcf().as_bytes())
                .unwrap();

            let header = vcf_reader.read_header().unwrap();

            let mut mapper = MapBasedMapper(map);
            let mut adapter = VCFToPopulationsAdapter::new(&header, None, &mut mapper).unwrap();

            for record in vcf_reader.records() {
                let record = record.unwrap();
                adapter.add_record(&record).unwrap();
            }

            let (population_name_to_idx, counts) = adapter.build();
            assert_eq!(population_name_to_idx.get("A").unwrap(), &0);
            assert_eq!(population_name_to_idx.get("B").unwrap(), &1);

            {
                let first_pop = &counts[0];
                let mut pop_iter = first_pop.iter();
                let first_site = pop_iter.next().unwrap();
                assert_eq!(first_site.counts(), &[5, 4]);
                assert_eq!(first_site.total_alleles(), 9);

                let second_site = pop_iter.next().unwrap();
                assert_eq!(second_site.counts(), &[4, 4]);
                assert_eq!(second_site.total_alleles(), 9);
            }

            {
                let second_pop = &counts[1];
                let mut pop_iter = second_pop.iter();
                let first_site = pop_iter.next().unwrap();
                assert_eq!(first_site.counts(), &[6, 3]);
                assert_eq!(first_site.total_alleles(), 9);

                let second_site = pop_iter.next().unwrap();
                assert_eq!(second_site.counts(), &[3, 4]);
                assert_eq!(second_site.total_alleles(), 9);
            }
        }
    }

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
    fn watterson_theta() {
        let sites = vec![
            vec![Some(0), Some(0), Some(1), Some(1), Some(1), Some(2)],
            vec![Some(0), Some(1), Some(1), Some(1), Some(2), None],
        ]
        .into_iter()
        .map(|site| {
            site.into_iter()
                .map(|sam| sam.map(AlleleID::from))
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();

        let allele_counts = MultiSiteCounts::from_tabular(sites);
        let theta = WattersonTheta::from_iter_sites(allele_counts.iter());

        let mut expected = 0f64;
        for site in allele_counts.iter() {
            let num_variants = site.counts.iter().filter(|c| **c > 0).count();
            let total_samples = site.counts.iter().filter(|c| **c > 0).sum::<Count>();

            if num_variants > 1 {
                let numerator = (num_variants - 1) as f64;
                let denominator = (1..total_samples).map(|i| 1f64 / i as f64).sum::<f64>();
                expected += numerator / denominator;
            }
        }

        assert!((theta.as_raw() - expected).abs() < f64::EPSILON);
    }

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

    #[test]
    fn watterson_theta_from_random_data() {
        use rand::prelude::*;

        let mut rng = StdRng::seed_from_u64(54321);
        for ploidy in [1, 2, 3, 4] {
            let freqs = [0.25, 0.5, 0.25]; // fixed allele freqs per site
            let mut sites = vec![];
            // make 10 random sites.
            // No missing data, etc..
            for _ in 0..10 {
                let site =
                    crate::testing::testdata::random_site_rng(10, ploidy, &freqs, None, &mut rng);
                sites.push(site);
            }
            // convert to our normal format
            let counts = crate::testing::testdata::single_pop_counts(&mut sites.iter());
            let theta = WattersonTheta::from_iter_sites(counts.iter());
            let theta_naive =
                crate::testing::naivecalculations::watterson_theta(&mut sites.iter_mut());
            if theta_naive.is_nan() {
                assert!(theta.as_raw().is_nan())
            } else {
                assert!((theta.as_raw() - theta_naive).abs() <= 1e-10)
            }
        }
    }

    #[test]
    fn watterson_theta_from_random_data_with_missing_data() {
        use rand::prelude::*;

        let mut rng = StdRng::seed_from_u64(54321);

        for ploidy in [1, 2, 4] {
            for rate in [0.01, 0.1, 0.5, 0.9] {
                let freqs = [0.25, 0.5, 0.25]; // fixed allele freqs per site
                let mut sites = vec![];
                // make 10 random sites.
                // No missing data, etc..
                for _ in 0..10 {
                    let site = crate::testing::testdata::random_site_rng(
                        10,
                        ploidy,
                        &freqs,
                        Some(crate::testing::testdata::RandomSiteOptions {
                            missing_data_rate: Some(rate),
                        }),
                        &mut rng,
                    );
                    sites.push(site);
                }
                // convert to our normal format
                let counts = crate::testing::testdata::single_pop_counts(&mut sites.iter());
                // get the calcs
                let theta = WattersonTheta::from_iter_sites(counts.iter());
                let theta_naive =
                    crate::testing::naivecalculations::watterson_theta(&mut sites.iter_mut());
                // compare
                if theta_naive.is_nan() {
                    assert!(theta.as_raw().is_nan());
                } else {
                    assert!(
                        (theta.as_raw() - theta_naive).abs() <= 1e-10,
                        "{theta:?} != {theta_naive}"
                    );
                }
            }
        }
    }
}
