use popgen::stats::StatRepresentation;
use popgen::stats::UnpolarisedSiteStat;

fn make_data() -> popgen::SampleAlleleCounts {
    let counts = vec![vec![Some(0.into()), None]];
    popgen::SampleAlleleCounts::try_from_tabular(counts).unwrap()
}

#[test]
fn diversity() {
    let counts = make_data();
    let diversity =
        popgen::stats::Diversity::try_from_iter_sites(counts.iter_sample_set(0).unwrap()).unwrap();
    assert_eq!(diversity.as_raw(), 0.)
}

#[test]
fn thetaw() {
    let counts = make_data();
    let thetaw =
        popgen::stats::WattersonsTheta::try_from_iter_sites(counts.iter_sample_set(0).unwrap())
            .unwrap();
    assert_eq!(thetaw.as_raw(), 0.)
}
