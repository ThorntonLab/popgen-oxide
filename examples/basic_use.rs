use popgen::stats::StatRepresentation; // Bring in ability for stats to express their value as a
                                       // low-level type via `as_raw`
use popgen::stats::UnpolarisedSiteStat; // Bring in the trait for Diversity to be able to do
                                        // calculations

fn main() {
    let data: Vec<Vec<Option<popgen::AlleleID>>> = vec![
        vec![Some(0.into()), Some(2.into())],
        vec![Some(0.into()), Some(1.into())],
    ];
    let counts = popgen::SampleAlleleCounts::try_from_tabular(data).unwrap();
    let diversity =
        popgen::stats::Diversity::try_from_iter_sites(counts.iter_sample_set(0).unwrap()).unwrap();
    assert_eq!(diversity.as_raw(), 2.0)
}
