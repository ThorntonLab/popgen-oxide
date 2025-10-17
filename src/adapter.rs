#[cfg(feature = "noodles")]
pub mod vcf {
    use crate::{AlleleID, PopgenResult};
    use noodles::vcf::variant::record::samples::keys::key;
    use noodles::vcf::variant::record::samples::series::Value;
    use noodles::vcf::variant::record::samples::Sample;
    use noodles::vcf::{Header, Record};

    pub use noodles::vcf as noodles_vcf;

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
                .transpose()?
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
                        genotypes.push(entry?.0.map(AlleleID::from))
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
}
