//! Basic integration tests

use std::num::NonZeroUsize;

use popgen::AlleleID;
use popgen_htslib::bcf_record_to_genotypes_iter_adapter;
use rust_htslib::bcf;
use rust_htslib::bcf::Read;

#[test]
fn test_basic_vcf_input_iter() {
    use std::io::Write;

    static VCF_FILE: &str = r#"##fileformat=VCFv4.6
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##contig=<ID=chr0>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	s0	s1	s2	s3	s4	s5	s6	s7	s8	s9	s10	s11	s12	s13	s14	s15	s16	s17
chr0	1	.	A	C	.	.	.	GT	0/0	0/1	0/1	0/0	0/1	0/0	0/1	0/0	0/0	0/0	0/0	1/1	1/0	0/1	0/0	./.	0/1	0/0
chr0	1	.	G	A	.	.	.	GT	0/0	./.	0/1	0/0	0/1	0/1	0/0	0/0	0/0	0/0	0/0	0/0	0/1	./.	0/1	0/1	0/0	0/0"#;
    let mut tfile = tempfile::NamedTempFile::new().unwrap();
    tfile.write_all(VCF_FILE.as_bytes()).unwrap();
    let (_, tfile_path) = tfile.into_parts();
    let mut bcf = bcf::Reader::from_path(tfile_path.as_os_str()).expect("Error opening file.");
    let mut counts =
        popgen::SampleAlleleCounts::of_empty_sample_sets(NonZeroUsize::new(1).unwrap());
    for record_result in bcf.records() {
        let record = record_result.unwrap();
        let allele_id_iter = popgen_htslib::bcf_record_to_genotypes_iter_adapter(&record)
            .unwrap()
            .flatten();
        counts.add_site(allele_id_iter).unwrap();
    }
    assert_eq!(
        counts
            .iter_sample_set(0)
            .unwrap()
            .filter(|c| c.total_alleles() == 36)
            .count(),
        2
    );
    assert_eq!(counts.num_sites(), 2);
    let num_non_missing = counts
        .iter_sample_set(0)
        .unwrap()
        .take(1)
        .map(|a| a.counts().iter().sum::<i64>())
        .collect::<Vec<_>>()[0];
    assert_eq!(num_non_missing, 2 * 18 - 2);
    let num_non_missing = counts
        .iter_sample_set(0)
        .unwrap()
        .skip(1)
        .map(|a| a.counts().iter().sum::<i64>())
        .collect::<Vec<_>>()[0];
    assert_eq!(num_non_missing, 2 * 18 - 4);
    for i in counts.iter_sample_set(0).unwrap() {
        assert_eq!(i.counts().len(), 2)
    }
}
#[test]

// The first half of samples are sample set 1.
// The second half are sample set 2.
fn test_basic_vcf_input_iter_to_multi_sample_set() {
    use std::io::Write;

    static VCF_FILE: &str = r#"##fileformat=VCFv4.6
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##contig=<ID=chr0>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	s0	s1	s2	s3	s4	s5	s6	s7	s8	s9	s10	s11	s12	s13	s14	s15	s16	s17
chr0	1	.	A	C	.	.	.	GT	0/0	0/1	0/1	0/0	0/1	0/0	0/1	0/0	0/0	0/0	0/0	1/1	1/0	0/1	0/0	./.	0/1	0/0
chr0	1	.	G	A	.	.	.	GT	0/0	./.	0/1	0/0	0/1	0/1	0/0	0/0	0/0	0/0	0/0	0/0	0/1	./.	0/1	0/1	0/0	0/0"#;
    let mut tfile = tempfile::NamedTempFile::new().unwrap();
    tfile.write_all(VCF_FILE.as_bytes()).unwrap();
    let (_, tfile_path) = tfile.into_parts();
    let mut bcf = bcf::Reader::from_path(tfile_path.as_os_str()).expect("Error opening file.");
    let mut counts =
        popgen::SampleAlleleCounts::of_empty_sample_sets(NonZeroUsize::new(2).unwrap());
    for record_result in bcf.records() {
        let record = record_result.unwrap();
        let mut allele_counts = vec![vec![]; 2];
        let mut num_alleles = [0_i64; 2];

        let genotype_iter = popgen_htslib::bcf_record_to_genotypes_iter_adapter(&record).unwrap();
        for (i, g) in genotype_iter.enumerate() {
            let index = if i < 9 { 0_usize } else { 1 };
            for allele in g {
                if let Some(id) = allele {
                    // The need for this implies that the newtype is a bit odd
                    let id = usize::from(id);
                    if id + 1 > allele_counts[index].len() {
                        allele_counts[index].resize(id + 1, 0);
                    }
                    allele_counts[index][id] += 1;
                }
                num_alleles[index] += 1;
            }
        }
        assert!(num_alleles.iter().all(|&i| i == 18));
        counts
            .extend_sample_sets_from_site_pred(|u| (&allele_counts[u], num_alleles[u]))
            .unwrap();
        assert!(counts
            .iter_sample_set(0)
            .unwrap()
            .map(|ac| ac.total_alleles())
            .all(|i| i == 18));
        assert!(counts
            .iter_sample_set(0)
            .unwrap()
            .all(|ac| ac.counts().len() == 2));
    }
}

#[test]
fn allele_correctness() {
    use std::io::Write;

    static VCF_FILE: &str = r#"##fileformat=VCFv4.6
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##contig=<ID=chr0>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	s0	s1	s2
chr0	1	.	A	C	.	.	.	GT	./.	0/0	1/1
chr0	2	.	A	C	.	.	.	GT	.|.	0|0	1|1"#;
    let mut tfile = tempfile::NamedTempFile::new().unwrap();
    tfile.write_all(VCF_FILE.as_bytes()).unwrap();
    let (_, tfile_path) = tfile.into_parts();
    let mut bcf = bcf::Reader::from_path(tfile_path.as_os_str()).expect("Error opening file.");

    let records = bcf
        .records()
        .map(|r| {
            bcf_record_to_genotypes_iter_adapter(&r.unwrap())
                .unwrap()
                .map(|gt| gt.collect::<Vec<_>>())
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();

    assert_eq!(
        records[0],
        vec![
            vec![None, None],
            vec![Some(AlleleID::from(0)), Some(AlleleID::from(0))],
            vec![Some(AlleleID::from(1)), Some(AlleleID::from(1))],
        ]
    );

    assert_eq!(
        records[1],
        vec![
            vec![None, None],
            vec![Some(AlleleID::from(0)), Some(AlleleID::from(0))],
            vec![Some(AlleleID::from(1)), Some(AlleleID::from(1))],
        ]
    );
}
