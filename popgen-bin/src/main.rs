use indicatif::{ProgressBar, ProgressIterator};
use popgen::adapter::record_to_genotypes_adapter;
use popgen::stats::{GlobalStatistic, WattersonTheta};
use popgen::MultiSiteCounts;
use std::error::Error;

fn main() -> Result<(), Box<dyn Error>> {
    let make_reader = || noodles::vcf::io::reader::Builder::default()
        // .build_from_path("ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf");
        .build_from_path("fix-ALL.chrY.phase3_integrated_v2b.20130502.genotypes.vcf");

    let num_records = {
        let mut reader = make_reader()?;
        reader.records().take(10000).count()
    };

    let mut reader = make_reader()?;

    let header = reader.read_header()?;
    let num_samples = header.sample_names().iter().count();

    let ploidy = 2;

    let allele_counts = MultiSiteCounts::from_tabular(reader.records()
        .take(10000)
        .progress_with(ProgressBar::new(num_records as u64))
        // drop Err records; we may want to revisit this
        .filter_map(Result::ok)
        .map(|rec| record_to_genotypes_adapter(&header, rec, num_samples, ploidy)));

    dbg!(WattersonTheta::from_iter_sites(allele_counts.iter()));
    // dbg!(allele_counts);

    Ok(())
}