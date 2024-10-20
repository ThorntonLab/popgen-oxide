use noodles::vcf::variant::record::samples::keys::key;
use noodles::vcf::variant::record::samples::series::Value;
use noodles::vcf::variant::record::samples::Sample as SampleTrait;
use popgen::{AlleleCounts, AlleleID};
use std::error::Error;
use indicatif::{ProgressBar, ProgressIterator};

fn main() -> Result<(), Box<dyn Error>> {
    let make_reader = || noodles::vcf::io::reader::Builder::default()
        .build_from_path("ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf");
    // .build_from_path("fix-ALL.chrY.phase3_integrated_v2b.20130502.genotypes.vcf");

    let num_records = {
        let mut reader = make_reader()?;
        reader.records().count()
    };

    let mut reader = make_reader()?;

    let header = reader.read_header()?;
    let num_samples = header.sample_names().iter().count();

    let ploidy = 2;

    let allele_counts = AlleleCounts::from_tabular(reader.records()
        .progress_with(ProgressBar::new(num_records as u64))
        // drop Err records; we may want to revisit this
        .filter_map(Result::ok)
        .map(|record| {
            let mut genotypes = Vec::with_capacity(ploidy * num_samples);
            record.samples().iter()
                .for_each(|sample| {
                    let fetched_field = match sample
                        // get the GT field
                        .get(&header, key::GENOTYPE)
                        .transpose()
                        .ok()
                        // bail if underlying IO fails
                        .expect("couldn't read GT") {
                        // return nothing if field missing
                        None => {
                            for _ in 0..ploidy {
                                genotypes.push(None);
                            }
                            return;
                        }
                        Some(fetched) => match fetched {
                            // return nothing if value missing
                            None => {
                                for _ in 0..ploidy {
                                    genotypes.push(None);
                                }
                                return;
                            }
                            // if everything checks out, proceed to the next match statement
                            Some(value) => value
                        }
                    };

                    match fetched_field {
                        Value::Genotype(genotype) => {
                            for entry in genotype.iter() {
                                genotypes.push(entry.expect("io error").0.map(AlleleID::from))
                            }
                        }
                        other => {
                            dbg!(other);
                            unreachable!()
                        }
                    };
                });

            genotypes.into_iter()
        }));

    dbg!(allele_counts);

    Ok(())
}