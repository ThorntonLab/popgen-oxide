#[cfg(feature = "noodles")]
pub mod vcf {
    use crate::{AlleleID, Count, PopgenResult, SampleAlleleCounts};
    use noodles::vcf::variant::record::samples::keys::key;
    use noodles::vcf::variant::record::samples::series::Value;
    use noodles::vcf::variant::record::samples::Sample;
    use noodles::vcf::variant::record::AlternateBases;
    use noodles::vcf::{Header, Record};
    use std::num::NonZeroI64;
    use std::ops::ControlFlow;

    pub fn record_to_genotypes_adapter(
        header: &Header,
        record: &Record,
        ploidy: usize,
    ) -> PopgenResult<Vec<Option<AlleleID>>> {
        let num_samples = header.sample_names().len();
        let mut genotypes = Vec::with_capacity(ploidy * num_samples);

        for sample in record.samples().iter() {
            let fetched_field = match sample
                // get the GT field
                .get(header, key::GENOTYPE)
                .transpose()
                .map_err(crate::PopgenError::NoodlesVCF)?
            {
                // return nothing if field missing
                None => {
                    for _ in 0..ploidy {
                        genotypes.push(None);
                    }
                    continue;
                }
                // return nothing if value missing
                Some(None) => {
                    for _ in 0..ploidy {
                        genotypes.push(None);
                    }
                    continue;
                }
                // if everything checks out, proceed to the next match statement
                Some(Some(value)) => value,
            };

            match fetched_field {
                Value::Genotype(genotype) => {
                    for entry in genotype.iter() {
                        genotypes.push(
                            entry
                                .map_err(crate::PopgenError::NoodlesVCF)?
                                .0
                                .map(AlleleID::from),
                        )
                    }
                }
                other => {
                    // panic because this has basically no reason to happen
                    dbg!(other);
                    panic!("parsed a genotype field and didn't get a genotype enum variant!");
                }
            };
        }
        Ok(genotypes)
    }

    /// `ploidy`, if not passed, will be inferred from the first record seen.
    pub struct VCFToSampleSetAdapter<'h> {
        header: &'h Header,
        ploidy: Option<NonZeroI64>,
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
        /// - `ploidy`: The ploidy in the data, or `None` to attempt to infer it from the first sample seen.
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
            ploidy: Option<NonZeroI64>,
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
                ploidy,
                sample_to_sample_set,
                sample_sets: SampleAlleleCounts::of_empty_sample_sets(num_sample_sets),
                // we'll resize if we ever get a record with more variants
                buf_counts: vec![0; num_sample_sets * 2],
                buf_num_samples: vec![0; num_sample_sets].into_boxed_slice(),
            })
        }

        pub fn add_record(&mut self, record: &Record) -> PopgenResult<()> {
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

            for (sample_i, sample) in record.samples().iter().enumerate() {
                let sample_set_id = self.sample_to_sample_set[sample_i];
                match sample
                    // get the GT field
                    .get(self.header, key::GENOTYPE)
                    .transpose()?
                    .flatten()
                {
                    // return nothing if field or value missing
                    None => {
                        let Some(ref ploidy) = self.ploidy else {
                            todo!("can't infer ploidy")
                        };

                        self.buf_num_samples[sample_set_id] += ploidy.get();
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

            self.sample_sets
                .extend_sample_sets_from_site(|sample_set_i| {
                    (
                        &self.buf_counts
                            [sample_set_i * num_sample_sets..(sample_set_i + 1) * num_sample_sets],
                        self.buf_num_samples[sample_set_i],
                    )
                })?;

            Ok(())
        }

        pub fn build(self) -> SampleAlleleCounts {
            self.sample_sets
        }
    }
}
