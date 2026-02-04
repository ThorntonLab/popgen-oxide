#[cfg(feature = "noodles")]
pub mod vcf {
    use crate::counts::MultiPopulationCounts;
    use crate::{AlleleID, PopgenResult};
    pub use noodles::vcf as noodles_vcf;
    use noodles::vcf::variant::record::samples::keys::key;
    use noodles::vcf::variant::record::samples::series::Value;
    use noodles::vcf::variant::record::samples::Sample;
    use noodles::vcf::variant::record::AlternateBases;
    use noodles::vcf::{Header, Record};
    use std::borrow::Cow;
    use std::collections::HashMap;
    use std::ops::ControlFlow;

    pub fn record_to_genotypes_adapter(
        header: &Header,
        record: Record,
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
    pub struct VCFToPopulationsAdapter<'h> {
        header: &'h Header,
        ploidy: Option<usize>,
        population_name_to_idx: HashMap<String, usize>,
        sample_to_population: Vec<usize>,
        populations: MultiPopulationCounts,
        // buffers for add_record
        buf_counts: Vec<i64>,
        buf_num_samples: Box<[usize]>,
    }

    pub trait WhichPopulation<'sample, 'pop, E> {
        fn which_population<'b>(&'b self, sample_name: &'sample str) -> Result<Cow<'pop, str>, E>
        where
            Self: 'pop,
            'b: 'pop;
    }

    impl<'sample, 'pop, E, T> WhichPopulation<'sample, 'pop, E> for T
    where
        T: Fn(&'sample str) -> Result<Cow<'pop, str>, E>,
    {
        fn which_population<'b>(&'b self, sample_name: &'sample str) -> Result<Cow<'pop, str>, E>
        where
            Self: 'pop,
            'b: 'pop,
        {
            self(sample_name)
        }
    }

    impl<'h> VCFToPopulationsAdapter<'h> {
        pub fn new<'map, 'sample, 'pop, W, E>(
            header: &'h Header,
            ploidy: Option<usize>,
            mapper: &'map W,
        ) -> Result<Self, E>
        where
            W: WhichPopulation<'sample, 'pop, E>,
            'h: 'sample,
            'map: 'pop,
        {
            let num_samples = header.sample_names().len();
            let mut sample_to_population = Vec::with_capacity(num_samples);
            let mut population_name_to_idx = HashMap::new();

            if let ControlFlow::Break(err) =
                header.sample_names().iter().try_for_each(|sample_name| {
                    sample_to_population.push({
                        let pop_name = match mapper.which_population(sample_name) {
                            Ok(name) => name,
                            Err(e) => return ControlFlow::Break(e),
                        };

                        if !population_name_to_idx.contains_key(&*pop_name) {
                            let population_id = population_name_to_idx.len();
                            population_name_to_idx.insert(pop_name.into_owned(), population_id);
                            population_id
                        } else {
                            *population_name_to_idx.get(&*pop_name).unwrap()
                        }
                    });

                    ControlFlow::Continue(())
                })
            {
                return Err(err);
            };

            let num_populations = population_name_to_idx.len();

            Ok(Self {
                header,
                ploidy,
                sample_to_population,
                population_name_to_idx,
                populations: MultiPopulationCounts::of_empty_populations(num_populations),
                // we'll resize if we ever get a record with more variants
                buf_counts: vec![0; num_populations * 2],
                buf_num_samples: vec![0; num_populations].into_boxed_slice(),
            })
        }

        pub fn add_record(&mut self, record: &Record) -> PopgenResult<()> {
            let num_populations = self.populations.num_populations();

            // let's assume that every stated allele is used
            let num_variants = 1 + record.alternate_bases().iter().count();

            let new_buf_counts_len = num_populations * num_variants;
            if new_buf_counts_len > self.buf_counts.len() {
                self.buf_counts.fill(0);
                self.buf_counts.resize(new_buf_counts_len, 0);
            } else {
                self.buf_counts.truncate(new_buf_counts_len);
                self.buf_counts.fill(0);
            }

            self.buf_num_samples.fill(0);

            for (sample_i, sample) in record.samples().iter().enumerate() {
                let population_id = self.sample_to_population[sample_i];
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

                        self.buf_num_samples[population_id] += *ploidy;
                    }
                    Some(Value::Genotype(genotype)) => {
                        for entry in genotype.iter() {
                            let (allele_id, _) = entry?;
                            self.buf_num_samples[population_id] += 1;
                            if let Some(allele_id) = allele_id {
                                self.buf_counts[population_id * num_populations + allele_id] += 1;
                            }
                        }
                    }
                    Some(_) => todo!("not a gt?"),
                };
            }

            self.populations
                .extend_populations_from_site(|population_i| {
                    (
                        &self.buf_counts
                            [population_i * num_populations..(population_i + 1) * num_populations],
                        self.buf_num_samples[population_i],
                    )
                })?;

            Ok(())
        }

        pub fn build(self) -> (HashMap<String, usize>, MultiPopulationCounts) {
            (self.population_name_to_idx, self.populations)
        }
    }
}
