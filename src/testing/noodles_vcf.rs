use crate::adapter::vcf::{record_to_genotypes_adapter, VCFToPopulationsAdapter, WhichPopulation};
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
    let mut writer = noodles::vcf::io::writer::Builder::default().build_from_writer(&mut buf);

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

fn counts_from_vcf(vcf_buf: &str, ploidy: usize) -> (Vec<Vec<Option<AlleleID>>>, MultiSiteCounts) {
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

    impl<'sample, 'pop> WhichPopulation<'sample, 'pop, ()> for MapBasedMapper {
        fn which_population<'b>(&'b self, sample_name: &'sample str) -> Result<Cow<'pop, str>, ()>
        where
            Self: 'pop,
            'b: 'pop,
        {
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

    let mapper = MapBasedMapper(map);
    let mut adapter = VCFToPopulationsAdapter::new(&header, None, &mapper).unwrap();

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

/// this test is mostly a check that this compiles
#[test]
fn load_vcf_multi_population_from_closure() {
    let mut vcf_reader = noodles::vcf::io::reader::Builder::default()
        .build_from_reader(make_vcf().as_bytes())
        .unwrap();

    let header = vcf_reader.read_header().unwrap();

    let mut adapter = VCFToPopulationsAdapter::new(&header, None, &|sample: &str| {
        let numeric_part = sample.split_at(1).1;
        let parsed = numeric_part.parse::<u8>().unwrap();
        let ret = match parsed % 2 {
            0 => Cow::Borrowed("A"),
            1 => Cow::Borrowed("B"),
            _ => unreachable!("mod 2 is either 0 or 1"),
        };

        Ok::<_, ()>(ret)
    }).unwrap();

    for record in vcf_reader.records() {
        let record = record.unwrap();
        adapter.add_record(&record).unwrap();
    }

    let (population_name_to_idx, counts) = adapter.build();
    assert_eq!(population_name_to_idx.get("A").unwrap(), &0);
    assert_eq!(population_name_to_idx.get("B").unwrap(), &1);
}
