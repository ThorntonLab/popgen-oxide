#[cfg(feature = "noodles")]
pub mod vcf {
    use crate::AlleleID;
    use noodles::vcf::variant::record::samples::keys::key;
    use noodles::vcf::variant::record::samples::series::Value;
    use noodles::vcf::variant::record::samples::Sample;
    use noodles::vcf::{Header, Record};

    pub use noodles::vcf as noodles_vcf;

    pub fn record_to_genotypes_adapter(
        header: &Header,
        record: Record,
        num_samples: usize,
        ploidy: usize,
    ) -> Vec<Option<AlleleID>> {
        let mut genotypes = Vec::with_capacity(ploidy * num_samples);
        record.samples().iter().for_each(|sample| {
            let fetched_field = match sample
                // get the GT field
                .get(header, key::GENOTYPE)
                .transpose()
                // bail if underlying IO fails
                .expect("couldn't read GT")
            {
                // return nothing if field missing
                None => {
                    for _ in 0..ploidy {
                        genotypes.push(None);
                    }
                    return;
                }
                // return nothing if value missing
                Some(None) => {
                    for _ in 0..ploidy {
                        genotypes.push(None);
                    }
                    return;
                }
                // if everything checks out, proceed to the next match statement
                Some(Some(value)) => value,
            };

            match fetched_field {
                Value::Genotype(genotype) => {
                    for entry in genotype.iter() {
                        genotypes.push(entry.expect("io error").0.map(AlleleID::from))
                    }
                }
                other => {
                    dbg!(other);
                    panic!("parsed a genotype field and didn't get a genotype enum variant!");
                }
            };
        });

        genotypes
    }
}
