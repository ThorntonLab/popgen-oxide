//! Adapter types for [`rust_htslib`]

use popgen::AlleleID;

/// Error type
#[non_exhaustive]
#[derive(Debug)]
pub enum Error {
    // NOTE: this is a bad name...
    /// Encapsulation of errors from [`rust_htslib`]
    RecordError(rust_htslib::errors::Error),
}

impl std::fmt::Display for Error {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::RecordError(e) => write!(f, "{e}"),
        }
    }
}

impl std::error::Error for Error {}

impl From<rust_htslib::errors::Error> for Error {
    fn from(value: rust_htslib::errors::Error) -> Self {
        Self::RecordError(value)
    }
}

/// Convert a BCF/VCF into allele counts
pub fn bcf_record_to_genotypes_adapter(
    record: &rust_htslib::bcf::Record,
) -> Result<Vec<Option<AlleleID>>, Error> {
    let mut site_counts_from_record = Vec::<Option<AlleleID>>::default();
    // NOTE: the error type is from std::num and we don't want that held in our
    // error type, so we need to map it to something else
    let sample_count = usize::try_from(record.sample_count()).unwrap();
    let gts = record.genotypes()?;
    for sample_index in 0..sample_count {
        for gt in gts.get(sample_index).iter() {
            match gt {
                rust_htslib::bcf::record::GenotypeAllele::Phased(index)
                | rust_htslib::bcf::record::GenotypeAllele::Unphased(index) => {
                    site_counts_from_record.push(Some((*index as usize).into()))
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
