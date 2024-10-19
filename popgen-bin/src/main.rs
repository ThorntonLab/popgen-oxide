use noodles::vcf::record::samples::Sample;
use noodles::vcf::variant::record::samples::keys::key;
use noodles::vcf::variant::record::samples::series::Value;
use noodles::vcf::variant::record::samples::Sample as SampleTrait;
use noodles::vcf::{Header, Record};
use popgen::AlleleID;
use std::error::Error;
use std::iter;

fn handle_sample<'a>(header: &'a Header, sample: Sample<'a>) -> Box<dyn Iterator<Item = Option<AlleleID>> + 'a> {
    match match sample
        // get the GT field
        .get(&header, key::GENOTYPE)
        .transpose()
        .ok()
        // bail if underlying IO fails
        .expect("couldn't read GT") {
        // return nothing if field missing
        None => return Box::new(iter::empty()),
        Some(fetched) => match fetched {
            // return nothing if value missing
            None => return Box::new(iter::empty()),
            // if everything checks out, proceed to the next match statement
            Some(value) => value
        }
    } {
        Value::Genotype(st) => Box::new(st.iter()
            .map(|res| res
                .expect("io error")
                .0.map(AlleleID::from))),
        other => {
            dbg!(other);
            panic!()
        }
    }
}

fn handle_record(header: &Header, record: Record) -> Box<dyn Iterator<Item = Option<AlleleID>>> {
    Box::from(record.samples().iter()
        .map(|sample| handle_sample(header, sample))
        .flatten()
    )
}

fn main() -> Result<(), Box<dyn Error>> {
    let mut reader = noodles::vcf::io::reader::Builder::default()
        // .build_from_path("ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf")?;
        .build_from_path("fix-ALL.chrY.phase3_integrated_v2b.20130502.genotypes.vcf")?;

    let header = reader.read_header()?;

    let count = reader.records()
        // drop Err records; we may want to revisit this
        .filter_map(Result::ok)
        .map(|record| handle_record(&header, record))
        .collect::<Vec<Box<dyn Iterator<Item = Option<AlleleID>>>>>();

    Ok(())
}