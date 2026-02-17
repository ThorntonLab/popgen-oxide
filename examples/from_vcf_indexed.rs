use noodles_core::Position;
use noodles_core::Region;
use popgen::adapter::vcf::record_to_genotypes_adapter;
use popgen::MultiSiteCounts;
use std::io::Write;

/*
Since this crate is agnostic to input data format, there is no special handling for e.g. indexed VCF.
Nonetheless, this example demonstrates that the adapter for VCF records also works for index-based queries.
*/

// these files are based on the following:
static VCF_FILE: &str = r#"##fileformat=VCFv4.5
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##contig=<ID=chr0>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	s0	s1	s2	s3	s4	s5	s6	s7
chr0	1	.	A	T	.	.	.	GT	/0	/1	/1	/0	/1	/0	/1	/0
chr0	2	.	G	C	.	.	.	GT	/0	/1	/0	/0	/1	/1	/0	/0
chr0	3	.	A	G	.	.	.	GT	/0	/1	/1	/1	/1	/1	/0	/1
chr0	4	.	C	G	.	.	.	GT	/0	/1	/1	/0	/1	/0	/0	/0"#;

fn main() {
    // Make a BGzipped file of the data shown above in VCF_FILE
    {
        let mut vcf_path = std::fs::File::create("simple_vcf.bgzf")
            .map(noodles::bgzf::Writer::new)
            .unwrap();
        vcf_path.write_all(VCF_FILE.as_bytes()).unwrap();
    }

    // Create an in-memory tabix-like index
    let index = noodles::vcf::fs::index("simple_vcf.bgzf").unwrap();
    let mut reader = noodles::vcf::io::indexed_reader::Builder::default()
        .set_index(index)
        .build_from_path("simple_vcf.bgzf")
        .unwrap();
    let header = reader.read_header().unwrap();
    let ploidy = 1;

    let region = Region::new(
        "chr0",
        Position::try_from(2).unwrap()..=Position::try_from(4).unwrap(),
    );

    let query = reader.query(&header, &region).unwrap();
    let alleles = query
        .map(Result::unwrap)
        .map(|rec| record_to_genotypes_adapter(&header, rec, ploidy).unwrap());
    let counts = MultiSiteCounts::from_tabular(alleles);
    counts.iter().for_each(|c| println!("{c:?}"));

    // clean up our file
    std::fs::remove_file("simple_vcf.bgzf").unwrap();
}
