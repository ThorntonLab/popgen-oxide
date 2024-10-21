use indicatif::{ProgressBar, ProgressIterator};
use popgen::adapter::record_to_genotypes_adapter;
use popgen::AlleleCounts;
use std::error::Error;

fn main() -> Result<(), Box<dyn Error>> {
    let make_reader = || noodles::vcf::io::reader::Builder::default()
        .build_from_path("ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf");
    // .build_from_path("fix-ALL.chrY.phase3_integrated_v2b.20130502.genotypes.vcf");

    let num_records = {
        let mut reader = make_reader()?;
        reader.records().count()
    };

    let mut reader = make_reader()?;

    let header = reader.read_header()?;
    let num_samples = header.sample_names().iter().count();

    let ploidy = 2;

    let allele_counts = AlleleCounts::from_tabular(reader.records()
        .progress_with(ProgressBar::new(num_records as u64))
        // drop Err records; we may want to revisit this
        .filter_map(Result::ok)
        .map(|rec| record_to_genotypes_adapter(&header, rec, num_samples, ploidy)));

    dbg!(allele_counts);

    Ok(())
}