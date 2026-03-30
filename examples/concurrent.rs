//! Example of concurrent computation enabled by the `SiteComposable` trait.

use noodles::vcf;
use popgen::adapter::vcf::record_to_genotypes_adapter;
use popgen::stats::{GlobalPi, GlobalStatistic, SiteComposable};
use popgen::MultiSiteCounts;
use rayon::iter::ParallelBridge;
use rayon::iter::ParallelIterator;
use std::io::{BufReader, Cursor};

static VCF_FILE: &str = r#"##fileformat=VCFv4.5
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##contig=<ID=chr0>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	s0	s1	s2	s3	s4	s5	s6	s7	s8	s9	s10	s11	s12	s13	s14	s15	s16	s17
chr0	1	.	A	C	.	.	.	GT	/0	/1	/1	/0	/1	/0	/1	/0	/0	/0	/0	/0	/0	/1	/0	/1	/1	/0
chr0	1	.	G	A	.	.	.	GT	/0	/1	/1	/0	/1	/1	/0	/0	/.	/.	/0	/0	/1	/1	/1	/1	/0	/."#;

fn main() {
    let mut reader = vcf::io::Reader::new(Cursor::new(VCF_FILE.as_bytes()));
    let header = reader.read_header().unwrap();

    let mut record = vcf::Record::default();
    let iter = std::iter::from_fn(move || match reader.read_record(&mut record).unwrap() {
        0 => None,
        _ => Some(record.clone()),
    });

    let all_alleles = iter
        .par_bridge()
        // this crate maps a record to allele IDs for you
        .map(|rec| record_to_genotypes_adapter(&header, &rec, 1).unwrap())
        .map(|genotypes| MultiSiteCounts::try_from_tabular(std::iter::once(genotypes)).unwrap())
        .map(|site| <GlobalPi as SiteComposable>::component_from(site.get(0).unwrap()))
        .fold(GlobalPi::default, |mut pi, c| {
            pi.try_add_component(c).unwrap();
            pi
        })
        // TODO: this is ugly
        .map(|pi| pi.as_raw())
        .reduce(f64::default, |a, b| a + b);

    dbg!(all_alleles);
}
