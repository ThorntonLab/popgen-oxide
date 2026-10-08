//! Adapter types for input via the [`noodles`] crate.

use noodles::vcf::variant::record::samples::keys::key;
use noodles::vcf::variant::record::samples::series::Value;
use noodles::vcf::variant::record::samples::Sample;
use noodles::vcf::variant::record::AlternateBases;
use noodles::vcf::{Header, Record};
use popgen::{AlleleID, Count, SampleAlleleCounts};
use std::ops::ControlFlow;

#[non_exhaustive]
#[derive(Debug)]
/// Error type
pub enum Error {
    /// Errors arising from the noodles crate when processing VCF records
    NoodlesVCF(std::io::Error),
    /// An input noodles [`Record`](Record) is badly formatted
    MalformedRecord,
    /// Contains [`popgen::PopgenError`]
    Popgen(popgen::PopgenError),
}

impl From<std::io::Error> for Error {
    fn from(e: std::io::Error) -> Self {
        Error::NoodlesVCF(e)
    }
}

impl From<popgen::PopgenError> for Error {
    fn from(e: popgen::PopgenError) -> Self {
        Error::Popgen(e)
    }
}

impl std::fmt::Display for Error {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Error::NoodlesVCF(e) => write!(f, "couldn't handle VCF: {}", e),
            Error::MalformedRecord => write!(f, "malformed VCF record"),
            Error::Popgen(e) => write!(f, "{e:?}"),
        }
    }
}

impl std::error::Error for Error {}

/// Process noodles input to a vector of [`Option`]al [`crate::AlleleID`]
///
/// # Parameters
///
/// * `header`: [`noodles::vcf::Header`]
/// * `record`: [`noodles::vcf::Record`]
/// * `ploidy`: the ploidy of the VCF records
pub fn record_to_alleles_adapter(
    header: &Header,
    record: &Record,
) -> Result<Vec<Option<AlleleID>>, Error> {
    let num_samples = header.sample_names().len();
    let mut num_genotypes_parsed = 0_usize;
    let mut genotypes = Vec::with_capacity(2 * num_samples);

    for sample in record.samples().iter() {
        let fetched_field = match sample
            // get the GT field
            .get(header, key::GENOTYPE)
            .transpose()
            .map_err(Error::NoodlesVCF)?
        {
            // return nothing if field missing
            None => {
                return Err(Error::MalformedRecord);
            }
            // return nothing if value missing
            Some(None) => {
                // This variant is not reachable.
                // Missing data gets handled below
                // in the match block.
                unreachable!("while parsing genotypes with noodles, the genotype field is present but the value is missing; this is a violation of the VCF spec");
            }
            // if everything checks out, proceed to the next match statement
            Some(Some(value)) => value,
        };

        match fetched_field {
            Value::Genotype(genotype) => {
                num_genotypes_parsed += 1;
                for entry in genotype.iter() {
                    // Here, the .0 is Option<usize>, and None implies missing data
                    genotypes.push(entry.map_err(Error::NoodlesVCF)?.0.map(AlleleID::from))
                }
            }
            other => {
                // panic because this has basically no reason to happen
                dbg!(other);
                panic!("parsed a genotype field and didn't get a genotype enum variant!");
            }
        };
    }
    if num_genotypes_parsed == num_samples {
        Ok(genotypes)
    } else {
        Err(Error::MalformedRecord)
    }
}

/// Builds [`popgen::SampleAlleleCounts`] from VCF
/// records for one or more sample sets.
pub struct VCFToSampleSetAdapter<'h> {
    header: &'h Header,
    sample_to_sample_set: Vec<usize>,
    sample_sets: SampleAlleleCounts,
    // buffers for add_record
    buf_counts: Vec<Count>,
    buf_num_samples: Box<[Count]>,
}

impl<'h> VCFToSampleSetAdapter<'h> {
    /// Build a new adapter.
    /// Requires:
    /// - `header`: A VCF header.
    /// - `num_sample_sets`: The number of sample sets.
    /// - `mapper`: An [`Fn`] from sample name (as `&str`) to a zero-based sample set ID.
    ///
    /// # Errors
    /// Any error from `mapper` will be propagated to the caller.
    ///
    /// # Panics
    /// If `mapper` produces a sample set ID greater than or equal to `num_sample_sets` (which is out-of-bounds in a zero-based ID system).
    pub fn new<'sample, M, E>(
        header: &'h Header,
        num_sample_sets: usize,
        mapper: M,
    ) -> Result<Self, E>
    where
        'h: 'sample,
        M: Fn(&'sample str) -> Result<usize, E>,
    {
        let num_samples = header.sample_names().len();
        let mut sample_to_sample_set = Vec::with_capacity(num_samples);

        if let ControlFlow::Break(err) =
                header.sample_names().iter().try_for_each(|sample_name| {
                    sample_to_sample_set.push(match mapper(sample_name) {
                        Ok(sample_set_id) if sample_set_id >= num_sample_sets => {
                            panic!("sample {sample_name} mapped to sample set ID {sample_set_id}, which is out of bounds for num_sample_sets {num_sample_sets}");
                        }
                        Ok(sample_set_id) => sample_set_id,
                        Err(e) => return ControlFlow::Break(e),
                    });

                    ControlFlow::Continue(())
                })
            {
                return Err(err);
            };

        Ok(Self {
            header,
            sample_to_sample_set,
            sample_sets: SampleAlleleCounts::of_empty_sample_sets(num_sample_sets),
            // we'll resize if we ever get a record with more variants
            buf_counts: vec![0; num_sample_sets * 2],
            buf_num_samples: vec![0; num_sample_sets].into_boxed_slice(),
        })
    }

    /// Process a [`noodles::vcf::Record`] into allele count data.
    pub fn add_record(&mut self, record: &Record) -> Result<(), Error> {
        let num_sample_sets = self.sample_sets.num_sample_sets();

        // let's assume that every stated allele is used
        let num_variants = 1 + record.alternate_bases().iter().count();

        let new_buf_counts_len = num_sample_sets * num_variants;
        if new_buf_counts_len > self.buf_counts.len() {
            self.buf_counts.fill(0);
            self.buf_counts.resize(new_buf_counts_len, 0);
        } else {
            self.buf_counts.truncate(new_buf_counts_len);
            self.buf_counts.fill(0);
        }

        self.buf_num_samples.fill(0);

        let mut samples_processed = 0_usize;
        for (sample_i, sample) in record.samples().iter().enumerate() {
            samples_processed += 1;
            let sample_set_id = self.sample_to_sample_set[sample_i];
            match sample
                // get the GT field
                .get(self.header, key::GENOTYPE)
                .transpose()?
                .flatten()
            {
                // return nothing if field or value missing
                None => {
                    return Err(Error::MalformedRecord);
                }
                Some(Value::Genotype(genotype)) => {
                    for entry in genotype.iter() {
                        let (allele_id, _) = entry?;
                        self.buf_num_samples[sample_set_id] += 1;
                        if let Some(allele_id) = allele_id {
                            self.buf_counts[sample_set_id * num_sample_sets + allele_id] += 1;
                        }
                    }
                }
                Some(_) => todo!("not a gt?"),
            };
        }
        if samples_processed != self.header.sample_names().len() {
            return Err(Error::MalformedRecord);
        }

        self.sample_sets
            .extend_sample_sets_from_site_pred(|sample_set_i| {
                (
                    &self.buf_counts
                        [sample_set_i * num_sample_sets..(sample_set_i + 1) * num_sample_sets],
                    self.buf_num_samples[sample_set_i],
                )
            })?;

        Ok(())
    }

    /// Consume `self`, returning [`crate::SampleAlleleCounts`]
    pub fn build(self) -> SampleAlleleCounts {
        self.sample_sets
    }
}
