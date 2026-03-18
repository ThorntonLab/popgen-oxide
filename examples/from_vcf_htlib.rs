use std::io::Write;

use rust_htslib::bcf;
use rust_htslib::bcf::Read;

static VCF_FILE: &str = r#"##fileformat=VCFv4.6
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##contig=<ID=chr0>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	s0	s1	s2	s3	s4	s5	s6	s7	s8	s9	s10	s11	s12	s13	s14	s15	s16	s17
chr0	1	.	A	C	.	.	.	GT	0/0	0/1	0/1	0/0	0/1	0/0	0/1	0/0	0/0	0/0	0/0	1/1	1/0	0/1	0/0	./.	0/1	0/0
chr0	1	.	G	A	.	.	.	GT	0/0	./.	0/1	0/0	0/1	0/1	0/0	0/0	0/0	0/0	0/0	0/0	0/1	./.	0/1	0/1	0/0	0/0"#;

fn main() {
    let mut vcf_out = std::fs::File::create_new("htslib_example.vcf").unwrap();
    vcf_out.write_all(VCF_FILE.as_bytes()).unwrap();
    let mut bcf = bcf::Reader::from_path("htslib_example.vcf").expect("Error opening file.");
    std::fs::remove_file("htslib_example.vcf").unwrap();
    let mut counts = popgen::MultiSiteCounts::default();
    let mut site_counts_from_record = Vec::<popgen::Count>::default();
    for record_result in bcf.records() {
        let record = record_result.unwrap();
        let num_alleles = record.alleles().len();
        site_counts_from_record.fill(0);
        site_counts_from_record.resize(num_alleles, 0);
        let gts = record.genotypes().expect("Error reading genotypes");
        // number of sample in the vcf
        let sample_count = usize::try_from(record.sample_count()).unwrap();
        let mut total_alleles = 0;
        for sample_index in 0..sample_count {
            // for each sample
            for gta in gts.get(sample_index).iter() {
                // Filter out missing genotypes
                if !matches!(gta, bcf::record::GenotypeAllele::PhasedMissing)
                    && !matches!(gta, bcf::record::GenotypeAllele::UnphasedMissing)
                {
                    site_counts_from_record[gta.index().unwrap() as usize] += 1;
                    total_alleles += 1;
                }
            }
        }
        counts
            .add_site_from_counts(site_counts_from_record.as_slice(), total_alleles)
            .unwrap();
    }
    println!("{counts:?}");
}
