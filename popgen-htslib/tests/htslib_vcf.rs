//! Basic integration tests

use rust_htslib::bcf;

#[test]
fn test_basic_vcf_input() {
    use std::io::Write;

    static VCF_FILE: &str = r#"##fileformat=VCFv4.6
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##contig=<ID=chr0>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	s0	s1	s2	s3	s4	s5	s6	s7	s8	s9	s10	s11	s12	s13	s14	s15	s16	s17
chr0	1	.	A	C	.	.	.	GT	0/0	0/1	0/1	0/0	0/1	0/0	0/1	0/0	0/0	0/0	0/0	1/1	1/0	0/1	0/0	./.	0/1	0/0
chr0	1	.	G	A	.	.	.	GT	0/0	./.	0/1	0/0	0/1	0/1	0/0	0/0	0/0	0/0	0/0	0/0	0/1	./.	0/1	0/1	0/0	0/0"#;
    let mut tfile = tempfile::tempfile().unwrap();
    tfile.write_all(VCF_FILE.as_bytes()).unwrap();
    let mut bcf = bcf::Reader::from_path(tfile).expect("Error opening file.");
    std::fs::remove_file("htslib_example.vcf").unwrap();
}
