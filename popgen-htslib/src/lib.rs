//! Adapter types for [`rust-htslib`]

use popgen::{AlleleID, PopgenResult};

pub fn bcf_record_to_genotypes_adapter(
    record: &rust_htslib::bcf::Record,
) -> PopgenResult<Vec<Option<AlleleID>>> {
    let mut site_counts_from_record = Vec::<Option<AlleleID>>::default();
    // NOTE: the error type is from std::num and we don't want that held in our
    // error type, so we need to map it to something else
    let sample_count = usize::try_from(record.sample_count()).unwrap();
    let gts = record.genotypes()?;
    for sample_index in 0..sample_count {
        for gt in gts.get(sample_index).iter() {
            match gt {
                rust_htslib::bcf::record::GenotypeAllele::Phased(_)
                | rust_htslib::bcf::record::GenotypeAllele::Unphased(_) => {
                    site_counts_from_record.push(Some(sample_index.into()))
                }
                rust_htslib::bcf::record::GenotypeAllele::PhasedMissing
                | rust_htslib::bcf::record::GenotypeAllele::UnphasedMissing => {
                    site_counts_from_record.push(None)
                }
            }
        }
    }
    Ok(site_counts_from_record)
}
