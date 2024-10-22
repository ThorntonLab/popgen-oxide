#[cfg(test)]
mod tests {
    use crate::noodles_vcf::variant::RecordBuf;
    use crate::{AlleleCounts, AlleleID};
    use noodles::vcf::header::record::value::map::{Contig, Format};
    use noodles::vcf::header::record::value::Map;
    use noodles::vcf::variant::io::Write;
    use noodles::vcf::variant::record::samples::keys::key;
    use noodles::vcf::variant::record_buf::samples::sample::value::Genotype;
    use noodles::vcf::variant::record_buf::samples::sample::Value;
    use noodles::vcf::variant::record_buf::samples::Keys;
    use noodles::vcf::variant::record_buf::{AlternateBases, Samples};

    use noodles::vcf::variant::record::samples::series::value::genotype::Phasing::Unphased;
    use noodles::vcf::variant::record_buf::samples::sample::value::genotype::Allele;
    use rand::seq::SliceRandom;
    use std::iter::repeat_n;

    // SiteVariant is to be a slice like ["A", "AG"] for a sample with these two genotypes
    // the appropriate IDs, number of samples, etc. will be calculated
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
            for gt_str in genotype.as_ref().iter()
                .filter_map(|gt_str| gt_str.as_ref())
                .map(|gt_str| gt_str.as_ref()) {
                if alleles_seen.iter().find(|allele| gt_str.eq(**allele)).is_none() {
                    alleles_seen.push(&*gt_str);
                }
            }
        }

        Some(RecordBuf::builder()
            .set_reference_sequence_name(seq_name)
            .set_variant_start(noodles_core::Position::MIN)
            .set_reference_bases(alleles_seen[0])
            .set_alternate_bases(AlternateBases::from(alleles_seen[1..].iter().map(|s| String::from(*s)).collect::<Vec<_>>()))
            .set_samples(Samples::new(
                Keys::from_iter(vec![String::from(key::GENOTYPE)]),
                // for each genotype and its count, create that many samples accordingly
                site_genotypes.iter().map(|(genotype, count)| vec![
                    Some(Value::from(Genotype::from_iter(
                        genotype.as_ref().iter()
                            .map(|sample_variant| Allele::new(
                                sample_variant.as_ref().map(|some| alleles_seen.iter()
                                    .enumerate()
                                    .find(|(_, variant)| **variant == some.as_ref())
                                    .unwrap().0
                                ),
                                Unphased  // TODO: make this configurable?
                            )),
                    )));
                    *count
                ])
                    .collect::<Vec<_>>()))
            .build()
        )
    }

    fn make_mock_vcf<Sites, SiteGenotypes, SiteVariant, Allele>(sites: Sites) -> Option<String>
    where
        Sites: AsRef<[SiteGenotypes]>,
        SiteGenotypes: AsRef<[(SiteVariant, usize)]>,
        SiteVariant: AsRef<[Option<Allele>]>,
        Allele: AsRef<str>,
    {
        let mut buf = Vec::new();
        let mut writer = noodles::vcf::io::writer::Builder::default()
            .build_from_writer(&mut buf);

        let seq_name = "chr0";
        // how many samples are needed to contain the genotypes of the sample with the most genotypes?
        let num_samples = sites.as_ref().iter()
            .map(|site_alleles| site_alleles.as_ref().iter()
                .map(|(_, count)| *count)
                .sum()
            ).max()?;

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

    #[test]
    fn load_raw() {
        let mut rng = rand::thread_rng();

        let sites = vec![
            {
                let mut site = repeat_n(0, 8)
                    .chain(repeat_n(1, 7))
                    .map(AlleleID::from)
                    .map(Option::from)
                    .chain(repeat_n(None, 4))
                    .collect::<Vec<_>>();
                site.shuffle(&mut rng);
                site
            },
            {
                let mut site = repeat_n(0, 341)
                    .chain(repeat_n(1, 69))
                    .chain(repeat_n(2, 926))
                    .map(AlleleID::from)
                    .map(Option::from)
                    .chain(repeat_n(None, 300))
                    .collect::<Vec<_>>();
                site.shuffle(&mut rng);
                site
            },
        ];

        let counts = AlleleCounts::from_tabular(sites);
        assert_eq!(counts.counts, vec![8, 7, 341, 69, 926]);
        assert_eq!(counts.count_starts, vec![0, 2]);
    }

    #[test]
    fn load_vcf() {
        let vcf_string = make_mock_vcf(vec![
            vec![
                (vec![Some("A")], 1),
                (vec![Some("T")], 24)
            ]
        ]).unwrap();

        dbg!(vcf_string);
    }
}
