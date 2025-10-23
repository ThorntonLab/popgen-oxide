use noodles_core::Position;
use noodles_core::Region;
use popgen::adapter::vcf::record_to_genotypes_adapter;
use popgen::MultiSiteCounts;
use std::io::Cursor;

/*
Since this crate is agnostic to input data format, there is no special handling for e.g. indexed VCF.
Nonetheless, this example demonstrates that the adapter for VCF records also works for index-based queries.
 */

// these files are based on the following:
#[allow(dead_code)]
static VCF_FILE: &str = r#"##fileformat=VCFv4.5
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##contig=<ID=chr0>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	s0	s1	s2	s3	s4	s5	s6	s7
chr0	1	.	A	T	.	.	.	GT	/0	/1	/1	/0	/1	/0	/1	/0
chr0	2	.	G	C	.	.	.	GT	/0	/1	/0	/0	/1	/1	/0	/0
chr0	3	.	A	G	.	.	.	GT	/0	/1	/1	/1	/1	/1	/0	/1
chr0	4	.	C	G	.	.	.	GT	/0	/1	/1	/0	/1	/0	/0	/0"#;

// bgzip simple.vcf
static VCF_BGZ: &[u8] = include_bytes!("simple.vcf.bgz");
// tabix simple.vcf.bgz
static VCF_TBI: &[u8] = include_bytes!("simple.vcf.bgz.tbi");

fn main() {
    let mut index_reader = noodles::tabix::io::Reader::new(VCF_TBI);

    let mut reader = noodles::vcf::io::indexed_reader::Builder::default()
        .set_index(index_reader.read_index().unwrap())
        // use Cursor to gain impl Seek for this example;
        // File is Seek so this would not be necessary in practice
        .build_from_reader(Cursor::new(VCF_BGZ))
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
}
