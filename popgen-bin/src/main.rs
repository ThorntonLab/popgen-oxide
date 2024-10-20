use itertools::Itertools;
use noodles::vcf::variant::record::samples::keys::key;
use noodles::vcf::variant::record::samples::series::Value;
use noodles::vcf::variant::record::samples::Sample as SampleTrait;
use popgen::{AlleleCounts, AlleleID};
use std::error::Error;


fn main() -> Result<(), Box<dyn Error>> {
    let mut reader = noodles::vcf::io::reader::Builder::default()
        // .build_from_path("ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf")?;
        .build_from_path("fix-ALL.chrY.phase3_integrated_v2b.20130502.genotypes.vcf")?;

    let header = reader.read_header()?;

    let ploidy = 2;

    let allele_counts = AlleleCounts::from_tabular(reader.records()
        // drop Err records; we may want to revisit this
        .filter_map(Result::ok)
        .map(|record| {
            record.samples().iter()
                .map(|sample| {
                    let fetched_field = match sample
                        // get the GT field
                        .get(&header, key::GENOTYPE)
                        .transpose()
                        .ok()
                        // bail if underlying IO fails
                        .expect("couldn't read GT") {
                        // return nothing if field missing
                        None => return vec![None; ploidy].into_iter(),
                        Some(fetched) => match fetched {
                            // return nothing if value missing
                            None => return vec![None; ploidy].into_iter(),
                            // if everything checks out, proceed to the next match statement
                            Some(value) => value
                        }
                    };

                    match fetched_field {
                        Value::Genotype(genotype) => {
                            genotype.iter()
                                .map(|res| res
                                    .expect("io error")
                                    .0.map(AlleleID::from))
                                .collect_vec()
                                .into_iter()
                        }
                        other => {
                            dbg!(other);
                            unreachable!()
                        }
                    }
                })
                .collect_vec()
                .into_iter()
        }));

    dbg!(allele_counts);

    Ok(())
}