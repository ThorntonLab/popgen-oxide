use popgen::stats::StatRepresentation;
use popgen::stats::UnpolarisedSiteStat;

// Test the case of a monomorphic site taking the form
// of one non-missing allele and the remainder are missing.

fn make_data_one_site() -> popgen::SampleAlleleCounts {
    let counts = vec![vec![Some(0.into()), None]];
    popgen::SampleAlleleCounts::try_from_tabular(counts).unwrap()
}

fn make_data_two_sites() -> popgen::SampleAlleleCounts {
    let counts = vec![vec![Some(0.into()), None]];
    let mut counts = popgen::SampleAlleleCounts::try_from_tabular(counts).unwrap();
    counts.add_site([Some(0.into()), Some(1.into())]).unwrap();
    counts
}

#[test]
fn diversity() {
    let counts = make_data_one_site();
    assert!(matches!(
        popgen::stats::Diversity::try_from_iter_sites(counts.iter_sample_set(0).unwrap()),
        Err(popgen::PopgenError::CalculationError)
    ));
    assert!(matches!(
        popgen::stats::Diversity::try_from_iter_sites(
            counts.iter_sample_set(0).unwrap().filter(|ac| ac
                .counts()
                .iter()
                .sum::<popgen::Count>()
                > 1),
        ),
        Err(popgen::PopgenError::EmptySiteCounts)
    ));
}

#[test]
fn diversity_two_sites() {
    let counts = make_data_two_sites();
    assert!(matches!(
        popgen::stats::Diversity::try_from_iter_sites(counts.iter_sample_set(0).unwrap()),
        Err(popgen::PopgenError::CalculationError)
    ));
    let div = popgen::stats::Diversity::try_from_iter_sites(
        counts
            .iter_sample_set(0)
            .unwrap()
            .filter(|ac| ac.counts().iter().sum::<popgen::Count>() > 1),
    )
    .unwrap();
    assert_eq!(div.as_raw(), 1.);
}

#[test]
fn thetaw() {
    let counts = make_data_one_site();
    let thetaw =
        popgen::stats::WattersonsTheta::try_from_iter_sites(counts.iter_sample_set(0).unwrap())
            .unwrap();
    assert_eq!(thetaw.as_raw(), 0.)
}

#[test]
fn tajd() {
    let counts = make_data_one_site();
    assert!(matches!(
        popgen::stats::TajimasD::try_from_iter_sites(counts.iter_sample_set(0).unwrap()),
        Err(popgen::PopgenError::CalculationError)
    ));
    assert!(matches!(
        popgen::stats::TajimasD::try_from_iter_sites(
            counts.iter_sample_set(0).unwrap().filter(|ac| ac
                .counts()
                .iter()
                .sum::<popgen::Count>()
                > 1),
        ),
        Err(popgen::PopgenError::EmptySiteCounts)
    ));
}
