use rust_htslib::bcf;
use rust_htslib::bcf::Read;

fn main() {
    let mut bcf = bcf::Reader::from_path("foo.vcf").expect("Error opening file.");
    let mut counts = popgen::MultiSiteCounts::default();
    let mut site_counts_from_record = Vec::<popgen::Count>::default();
    for record_result in bcf.records() {
        let record = record_result.unwrap();
        let num_alleles = record.alleles().len();
        site_counts_from_record.resize(num_alleles, 0);
        let gts = record.genotypes().expect("Error reading genotypes");
        // number of sample in the vcf
        let sample_count = usize::try_from(record.sample_count()).unwrap();
        let mut total_alleles = 0;
        for sample_index in 0..sample_count {
            // for each sample
            for gta in gts.get(sample_index).iter() {
                site_counts_from_record[gta.index().unwrap() as usize] += 1;
                total_alleles += 1; // FIXME: this is wrong...
            }
        }
        counts
            .add_site_from_counts(site_counts_from_record.as_slice(), total_alleles)
            .unwrap();
    }
}
