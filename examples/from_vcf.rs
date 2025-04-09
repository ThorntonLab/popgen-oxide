use noodles::vcf;
use popgen::{adapter::vcf::record_to_genotypes_adapter, MultiSiteCounts};

static VCF_FILE: &str = r#"##fileformat=VCFv4.5
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##contig=<ID=chr0>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	s0	s1	s2	s3	s4	s5	s6	s7	s8	s9	s10	s11	s12	s13	s14	s15	s16	s17
chr0	1	.	A	C	.	.	.	GT	/0	/1	/1	/0	/1	/0	/1	/0	/0	/0	/0	/0	/0	/1	/0	/1	/1	/0
chr0	1	.	G	A	.	.	.	GT	/0	/1	/1	/0	/1	/1	/0	/0	/.	/.	/0	/0	/1	/1	/1	/1	/0	/."#;

fn main() {
    let mut reader = vcf::io::reader::Builder::default()
        // or read something yourself. of course
        .build_from_reader(VCF_FILE.as_bytes())
        .unwrap();

    let ploidy = 1;

    let header = reader.read_header().unwrap();
    let num_samples = header.sample_names().len();
    let all_alleles = reader
        .records()
        // ignore IO errors
        .map(Result::unwrap)
        // this crate maps a record to allele IDs for you
        .map(|rec| record_to_genotypes_adapter(&header, rec, num_samples, ploidy).unwrap());

    // this constructor is iterator-based
    let counts = MultiSiteCounts::from_tabular(all_alleles);
    counts.iter().for_each(|c| println!("{c:?}"));
}
