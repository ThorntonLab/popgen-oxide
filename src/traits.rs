/// Types implementing this trait can be reduced by consuming two values to produce a new one.
///
/// For count types like [`crate::counts::SampleAlleleCounts`], this means concatenating sites from the two inputs.
/// For statistic types like [`crate::stats::Diversity`], this means composing the results from different data, as if the statistic was computed on the concatenated data originally.
pub trait TryReduce {
    type Error: std::error::Error;

    fn try_reduce(self, other: Self) -> Result<Self, Self::Error>
    where
        Self: Sized;
}

#[test]
fn sample_allele_counts_reduce() {
    use crate::AlleleID;
    use crate::Count;
    use crate::SampleAlleleCounts;

    let counts1 = SampleAlleleCounts::try_from_tabular(
        vec![[0, 0, 0], [0, 1, 1], [0, 1, 2]]
            .into_iter()
            .map(|site| site.into_iter().map(|a| Some(AlleleID::from(a)))),
    )
    .unwrap();

    let counts2 = SampleAlleleCounts::try_from_tabular(
        vec![[1, 1, 1], [2, 2, 2], [0, 0, 2]]
            .into_iter()
            .map(|site| site.into_iter().map(|a| Some(AlleleID::from(a)))),
    )
    .unwrap();

    let combined = counts1.try_reduce(counts2).unwrap();
    assert_eq!(combined.iter().count(), 6);
    let expect: Vec<(&[Count], i32)> = vec![
        (&[3], 3),
        (&[1, 2], 3),
        (&[1, 1, 1], 3),
        (&[0, 3], 3),
        (&[0, 0, 3], 3),
        (&[2, 0, 1], 3),
    ];

    for (e, a) in expect.into_iter().zip(combined.iter()) {
        assert_eq!(e.0, a.counts());
        assert_eq!(e.1, a.total_alleles());
    }
}
