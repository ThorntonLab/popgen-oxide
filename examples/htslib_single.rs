use rust_htslib::bcf::record::GenotypeAllele;
use rust_htslib::bcf::Read;
use popgen::stats::{GlobalPi, GlobalStatistic};

fn main() {
    let mut bcf = rust_htslib::bcf::IndexedReader::from_path(
        "ALL.chr18.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz",
    )
    .expect("Error opening file.");

    let mut counts = popgen::MultiSiteCounts::default();
    let mut site_counts_from_record = Vec::<popgen::Count>::default();

    for record_result in bcf.records() {
        let record = record_result.unwrap();
        let num_alleles = record.alleles().len();

        if num_alleles > site_counts_from_record.len() {
            site_counts_from_record.fill(0);
            site_counts_from_record.resize(num_alleles, 0);
        } else {
            site_counts_from_record.truncate(num_alleles);
            site_counts_from_record.fill(0);
        }

        let gts = record.genotypes().expect("Error reading genotypes");

        // number of sample in the vcf
        let sample_count = usize::try_from(record.sample_count()).unwrap();
        let mut total_alleles = 0;
        for sample_index in 0..sample_count {
            // for each sample
            for gta in gts.get(sample_index).iter() {
                // Filter out missing genotypes
                if let GenotypeAllele::Phased(ind) | GenotypeAllele::Unphased(ind) = gta {
                    site_counts_from_record[*ind as usize] += 1;
                }
                total_alleles += 1;
            }
        }
        counts
            .add_site_from_counts(site_counts_from_record.as_slice(), total_alleles)
            .unwrap();
    }

    dbg!(GlobalPi::try_from_iter_sites(counts.iter()));
}
