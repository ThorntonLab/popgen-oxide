#[cfg(test)]
pub fn is_close(a: f64, b: f64) -> bool {
    assert!(a.is_normal());
    assert!(b.is_normal());

    (a - b).abs() < f64::EPSILON
}

#[cfg(test)]
mod model {
    use popgen::{AlleleID, MultiSiteCounts};
    use std::collections::HashMap;

    pub type Allele = Option<Box<[u8]>>;
    pub type IndividualGenotype = Vec<Allele>;
    pub type Site = Vec<IndividualGenotype>;
    #[derive(Clone, Debug, Default)]
    pub struct GenomeCollection {
        data: Vec<Site>,
    }

    impl GenomeCollection {
        pub fn new(data: Vec<Site>) -> Self {
            Self { data }
        }

        pub fn from_iter_sites<SiteIt, IndividualIt, GenotypeIt>(iter: SiteIt) -> Self
        where
            SiteIt: Iterator<Item = IndividualIt>,
            IndividualIt: Iterator<Item = GenotypeIt>,
            GenotypeIt: Iterator<Item = Allele>,
        {
            Self {
                data: iter
                    .map(|indiv| indiv.map(|gt| gt.collect()).collect())
                    .collect(),
            }
        }

        pub fn sites(&self) -> impl Iterator<Item = &Site> {
            self.data.iter()
        }

        pub fn samples(&self) -> SamplesIter {
            SamplesIter {
                current_sample: 0,
                num_samples: self.data[0].len(),
                inner: self,
            }
        }
    }

    #[derive(Debug)]
    pub struct SamplesIter<'a> {
        current_sample: usize,
        num_samples: usize,
        inner: &'a GenomeCollection,
    }

    impl<'a> Iterator for SamplesIter<'a> {
        type Item = SampleSitesIter<'a>;

        fn next(&mut self) -> Option<Self::Item> {
            if self.current_sample >= self.num_samples {
                None
            } else {
                let ret = Some(SampleSitesIter {
                    current_sample: self.current_sample,
                    current_site: 0,
                    num_sites: self.inner.data.len(),
                    inner: self.inner,
                });
                self.current_sample += 1;
                ret
            }
        }
    }

    pub struct SampleSitesIter<'a> {
        current_sample: usize,
        current_site: usize,
        num_sites: usize,
        inner: &'a GenomeCollection,
    }

    impl<'a> Iterator for SampleSitesIter<'a> {
        type Item = &'a IndividualGenotype;

        fn next(&mut self) -> Option<Self::Item> {
            if self.current_site >= self.num_sites {
                None
            } else {
                let ret = Some(&self.inner.data[self.current_site][self.current_sample]);
                self.current_site += 1;
                ret
            }
        }
    }

    impl From<GenomeCollection> for MultiSiteCounts {
        fn from(val: GenomeCollection) -> Self {
            // convert sites to IDs; suppose the variants appearing first get lower IDs
            MultiSiteCounts::from_tabular(val.data.iter().map(|site| {
                let mut ret = Vec::with_capacity(site.iter().map(|indiv| indiv.len()).sum());
                let mut counts = HashMap::new();
                let mut next_id = 0usize;
                for individual in site {
                    for genotype in individual {
                        ret.push(genotype.as_ref().map(|actual| {
                            let assigned_id = counts.entry(&**actual).or_insert(next_id);
                            if assigned_id == &next_id {
                                next_id += 1;
                            }

                            AlleleID::from(*assigned_id)
                        }));
                    }
                }
                ret.into_iter()
            }))
        }
    }
}

#[cfg(test)]
mod naive_conversions {
    #[cfg(feature = "noodles")]
    mod noodles {
        use std::io::BufRead;

        use noodles::vcf::{
            variant::record::samples::{keys::key, series::Value, Sample},
            variant::record::AlternateBases,
            Header, Record,
        };
        use popgen::{adapter::vcf::record_to_genotypes_adapter, MultiSiteCounts, PopgenResult};

        use crate::model::{GenomeCollection, Site};

        impl GenomeCollection {
            fn vcf_record_to_site(header: &Header, record: &Record, ploidy: usize) -> Site {
                let mut ret = Vec::with_capacity(ploidy * record.samples().iter().count());

                for sample in record.samples().iter() {
                    let Some(gt_field) = sample
                        .get(header, key::GENOTYPE)
                        .transpose()
                        .unwrap()
                        .flatten()
                    else {
                        for _ in 0..ploidy {
                            ret.push(vec![]);
                        }
                        continue;
                    };

                    let Value::Genotype(genotype) = gt_field else {
                        panic!("genotype field is not a genotype-type value")
                    };

                    let individual = genotype
                        .iter()
                        .map(|entry| {
                            entry.unwrap().0.map(|id| match id {
                                0 => Box::<[u8]>::from(record.reference_bases().as_bytes()),
                                other => {
                                    let bases = record.alternate_bases();
                                    let nth = bases.iter().nth(other - 1);
                                    Box::<[u8]>::from(nth.unwrap().unwrap().as_bytes())
                                }
                            })
                        })
                        .collect::<Vec<_>>();

                    ret.push(individual);
                }

                ret
            }

            pub fn from_vcf<R: BufRead>(
                mut reader: noodles::vcf::io::Reader<R>,
                ploidy: usize,
            ) -> Self {
                let header = reader.read_header().unwrap();

                Self::new(
                    reader
                        .records()
                        .map(Result::unwrap)
                        .map(|rec| Self::vcf_record_to_site(&header, &rec, ploidy))
                        .collect::<Vec<Site>>(),
                )
            }
        }

        #[test]
        fn naive_vcf() {
            const VCF: &str = r#"##fileformat=VCFv4.5
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##contig=<ID=chr0>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	s0	s1	s2	s3	s4	s5	s6	s7	s8	s9	s10	s11	s12	s13	s14	s15	s16	s17
chr0	1	.	G	CT	.	.	.	GT	/0	/1	/1	/0	/1	/0	/1	/0	/0	/0	/1	/0	/0	/1	/0	/1	/1	/0
chr0	1	.	CG	A,T	.	.	.	GT	/0	/1	/1	/0	/1	/1	/2	/0	/.	/.	/0	/2	/1	/1	/1	/1	/0	/."#;

            let ploidy = 1;

            let from_naive: MultiSiteCounts = {
                let reader = noodles::vcf::io::reader::Builder::default()
                    .build_from_reader(VCF.as_bytes())
                    .unwrap();

                GenomeCollection::from_vcf(reader, ploidy)
            }
            .into();

            let counts = {
                let mut reader = noodles::vcf::io::reader::Builder::default()
                    .build_from_reader(VCF.as_bytes())
                    .unwrap();

                let header = reader.read_header().unwrap();
                let num_samples = header.sample_names().iter().count();

                let all_alleles = reader
                    .records()
                    .map(Result::unwrap)
                    .map(|rec| record_to_genotypes_adapter(&header, rec, num_samples, ploidy))
                    .collect::<PopgenResult<Vec<_>>>()
                    .unwrap();
                MultiSiteCounts::from_tabular(all_alleles)
            };

            assert!(from_naive.iter().zip(counts.iter()).all(|(a, b)| a == b));
        }
    }
}

#[cfg(test)]
mod naive_stats {
    use std::collections::HashSet;

    use itertools::repeat_n;
    use popgen::{
        stats::{GlobalPi, GlobalStatistic, WattersonTheta},
        MultiSiteCounts,
    };
    use rand::{rng, rngs::ThreadRng, seq::SliceRandom};
    use triangle_matrix::{
        SimpleLowerTri, SymmetricUpperTri, SymmetricUpperTriMut, Triangle, TriangleMut,
    };

    use crate::{
        is_close,
        model::{GenomeCollection, IndividualGenotype, Site},
    };

    trait NaiveGlobalStatistic: Default {
        fn as_raw(&self) -> f64;

        fn from_iter_sites<'a>(sites: impl Iterator<Item = &'a Site>) -> Self;

        fn from_iter_samples<'a>(
            samples: impl Iterator<Item = impl Iterator<Item = &'a IndividualGenotype>>,
        ) -> Self;
    }

    struct TriVec<T>(usize, Vec<T>);

    impl<T> Triangle<T> for TriVec<T> {
        type Inner = Vec<T>;

        fn n(&self) -> usize {
            self.0
        }

        fn inner(&self) -> &Self::Inner {
            &self.1
        }
    }

    impl<T> TriangleMut<T> for TriVec<T> {
        fn inner_mut(&mut self) -> &mut Self::Inner {
            &mut self.1
        }
    }

    #[derive(Debug, Default)]
    struct NaiveGlobalPi(f64);

    impl NaiveGlobalStatistic for NaiveGlobalPi {
        fn as_raw(&self) -> f64 {
            self.0
        }

        fn from_iter_sites<'a>(sites: impl Iterator<Item = &'a Site>) -> Self {
            let mut ret = Self::default();
            for site in sites {
                let total_alleles = site.iter().map(|indiv| indiv.len()).sum::<usize>();

                let mut mat = TriVec(
                    total_alleles,
                    Vec::from_iter(repeat_n(
                        None::<i32>,
                        total_alleles * (total_alleles + 1) / 2,
                    )),
                );

                for (ind_a, allele_a) in site.iter().flat_map(|indiv| indiv.iter()).enumerate() {
                    for (ind_b, allele_b) in site
                        .iter()
                        .flat_map(|indiv| indiv.iter())
                        .enumerate()
                        .skip(ind_a + 1)
                    {
                        *SymmetricUpperTriMut::get_element_mut(&mut mat, ind_a, ind_b) =
                            match *allele_a {
                                None => None,
                                Some(ref allele_id_a) => match allele_b {
                                    None => None,
                                    Some(allele_id_b) => match allele_id_a.eq(allele_id_b) {
                                        false => Some(1),
                                        true => Some(0),
                                    },
                                },
                            }
                    }
                }

                let make_iter = || {
                    mat.iter_triangle_indices()
                        .filter_map(|(i, j)| *SymmetricUpperTri::get_element(&mat, i, j))
                };
                ret.0 += make_iter().sum::<i32>() as f64 / make_iter().count() as f64;
            }
            ret
        }

        fn from_iter_samples<'a>(
            samples: impl Iterator<Item = impl Iterator<Item = &'a IndividualGenotype>>,
        ) -> Self {
            let mut total_differences = 0f64;
            let mut samples_collect = Vec::<Vec<IndividualGenotype>>::new();

            for s in samples {
                let new_sample = s.cloned().collect::<Vec<_>>();
                let total_alleles = new_sample.iter().map(|indiv| indiv.len()).sum::<usize>();

                let mut mat = TriVec(
                    total_alleles,
                    Vec::from_iter(repeat_n(
                        None::<i32>,
                        total_alleles * (total_alleles + 1) / 2,
                    )),
                );

                for (ind_a, allele_a) in
                    new_sample.iter().flat_map(|indiv| indiv.iter()).enumerate()
                {
                    for (ind_b, allele_b) in new_sample
                        .iter()
                        .flat_map(|indiv| indiv.iter())
                        .enumerate()
                        .skip(ind_a + 1)
                    {
                        *SymmetricUpperTriMut::get_element_mut(&mut mat, ind_a, ind_b) =
                            match *allele_a {
                                None => None,
                                Some(ref allele_id_a) => match allele_b {
                                    None => None,
                                    Some(allele_id_b) => match allele_id_a.eq(allele_id_b) {
                                        false => Some(1),
                                        true => Some(0),
                                    },
                                },
                            }
                    }
                }

                let make_iter = || {
                    mat.iter_triangle_indices()
                        .filter_map(|(i, j)| *SymmetricUpperTri::get_element(&mat, i, j))
                };
                total_differences += make_iter().sum::<i32>() as f64 / make_iter().count() as f64;

                samples_collect.push(new_sample);
            }

            Self(total_differences)
        }
    }

    fn box_allele<const N: usize>(allele: Option<&'static [u8; N]>) -> Option<Box<[u8]>> {
        allele.map(|a| Box::from(&a[..]))
    }

    fn shuffled_alleles<I: Clone>(
        ids: impl Iterator<Item = (I, usize)>,
        rng: &mut ThreadRng,
    ) -> Vec<I> {
        let mut site = vec![];
        ids.for_each(|(id, count)| {
            site.extend(repeat_n(id, count));
        });

        site.shuffle(rng);
        site
    }

    #[test]
    fn global_pi_against_naive_from_sites() {
        let mut rng = rng();
        let sites = vec![
            vec![shuffled_alleles(
                [(Some(b"AGA"), 35), (Some(b"GAG"), 6), (None, 3)]
                    .into_iter()
                    .map(|(a, b)| (box_allele(a), b))
                    .collect::<Vec<_>>()
                    .into_iter(),
                &mut rng,
            )],
            vec![shuffled_alleles(
                [(Some(b"GGT"), 2), (Some(b"AGC"), 14), (Some(b"ACC"), 155)]
                    .into_iter()
                    .map(|(a, b)| (box_allele(a), b))
                    .collect::<Vec<_>>()
                    .into_iter(),
                &mut rng,
            )],
        ];

        let collection = GenomeCollection::new(sites);
        let allele_counts: MultiSiteCounts = collection.clone().into();
        assert!(is_close(
            NaiveGlobalPi::from_iter_sites(collection.sites()).as_raw(),
            GlobalPi::from_iter_sites(allele_counts.iter()).as_raw()
        ));
    }

    #[test]
    fn global_pi_against_naive_against_samples() {
        let mut rng = rng();

        let sites = vec![
            vec![shuffled_alleles(
                [(Some(b"AGA"), 35), (Some(b"GAG"), 6), (None, 3)]
                    .into_iter()
                    .map(|(a, b)| (box_allele(a), b))
                    .collect::<Vec<_>>()
                    .into_iter(),
                &mut rng,
            )],
            vec![shuffled_alleles(
                [(Some(b"GGT"), 2), (Some(b"AGC"), 14), (Some(b"ACC"), 155)]
                    .into_iter()
                    .map(|(a, b)| (box_allele(a), b))
                    .collect::<Vec<_>>()
                    .into_iter(),
                &mut rng,
            )],
        ];

        let collection = GenomeCollection::new(sites);
        let allele_counts: MultiSiteCounts = collection.clone().into();
        assert!(is_close(
            NaiveGlobalPi::from_iter_samples(collection.samples()).as_raw(),
            GlobalPi::from_iter_sites(allele_counts.iter()).as_raw()
        ));
    }

    #[derive(Debug, Default)]
    struct NaiveWattersonTheta(f64);

    impl NaiveGlobalStatistic for NaiveWattersonTheta {
        fn as_raw(&self) -> f64 {
            self.0
        }

        fn from_iter_sites<'a>(sites: impl Iterator<Item = &'a Site>) -> Self {
            let mut ret = Self::default();

            for site in sites {
                let num_variants = site
                    .iter()
                    .flat_map(|indiv| indiv.iter())
                    .filter_map(|o| o.as_ref())
                    .collect::<HashSet<_>>()
                    .len();
                let total_samples = site
                    .iter()
                    .flat_map(|indiv| indiv.iter())
                    .filter_map(|o| o.as_ref())
                    .count();

                if num_variants > 1 {
                    let numerator = (num_variants - 1) as f64;
                    let denominator = (1..total_samples).map(|i| 1f64 / i as f64).sum::<f64>();
                    ret.0 += numerator / denominator;
                }
            }
            ret
        }

        fn from_iter_samples<'a>(
            samples: impl Iterator<Item = impl Iterator<Item = &'a IndividualGenotype>>,
        ) -> Self {
            let mut samples = samples.map(|site| site.collect::<Vec<_>>()).peekable();
            let Some(first_sample) = samples.peek() else {
                return Self::default();
            };

            let num_sites = first_sample.len();
            let mut state_appearance_per_site = vec![HashSet::new(); num_sites];
            let mut sample_count_per_site = vec![0; num_sites];

            for sample in samples {
                for (site_i, site) in sample.iter().enumerate() {
                    for allele in site.iter().cloned() {
                        if let Some(real) = allele {
                            sample_count_per_site[site_i] += 1;
                            state_appearance_per_site[site_i].insert(real);
                        }
                    }
                }
            }

            let mut ret = Self::default();
            for (appearance, count) in state_appearance_per_site.iter().zip(sample_count_per_site) {
                let s_i = appearance.len() - 1;

                let mut denom = 0f64;
                for i in 1..count {
                    denom += 1f64 / i as f64;
                }

                ret.0 += s_i as f64 / denom;
            }

            ret
        }
    }

    #[test]
    fn watterson_theta_against_naive_from_sites() {
        let mut rng = rng();

        let sites = vec![
            vec![shuffled_alleles(
                [(Some(b"A"), 2), (Some(b"T"), 3), (Some(b"G"), 1)]
                    .into_iter()
                    .map(|(a, b)| (box_allele(a), b))
                    .collect::<Vec<_>>()
                    .into_iter(),
                &mut rng,
            )],
            vec![shuffled_alleles(
                [(Some(b"T"), 1), (Some(b"C"), 3), (Some(b"G"), 1), (None, 1)]
                    .into_iter()
                    .map(|(a, b)| (box_allele(a), b))
                    .collect::<Vec<_>>()
                    .into_iter(),
                &mut rng,
            )],
        ];
        let collection = GenomeCollection::new(sites);
        let allele_counts: MultiSiteCounts = collection.clone().into();

        let theta_naive = NaiveWattersonTheta::from_iter_sites(collection.sites());
        let theta = WattersonTheta::from_iter_sites(allele_counts.iter());

        assert!(is_close(theta_naive.as_raw(), theta.as_raw()));
    }

    #[test]
    fn watterson_theta_against_naive_from_samples() {
        let mut rng = rng();

        let sites = vec![
            vec![shuffled_alleles(
                [(Some(b"A"), 2), (Some(b"T"), 3), (Some(b"G"), 1)]
                    .into_iter()
                    .map(|(a, b)| (box_allele(a), b))
                    .collect::<Vec<_>>()
                    .into_iter(),
                &mut rng,
            )],
            vec![shuffled_alleles(
                [(Some(b"T"), 1), (Some(b"C"), 3), (Some(b"G"), 1), (None, 1)]
                    .into_iter()
                    .map(|(a, b)| (box_allele(a), b))
                    .collect::<Vec<_>>()
                    .into_iter(),
                &mut rng,
            )],
        ];
        let collection = GenomeCollection::new(sites);
        let allele_counts: MultiSiteCounts = collection.clone().into();

        let theta_naive = NaiveWattersonTheta::from_iter_samples(collection.samples());
        let theta = WattersonTheta::from_iter_sites(allele_counts.iter());
        dbg!((&theta_naive, &theta));

        assert!(is_close(theta_naive.as_raw(), theta.as_raw()));
    }
}
