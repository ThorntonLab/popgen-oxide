#[test]
fn test_iteration_over_empty() {
    let counts = SampleAlleleCounts::default();
    assert_eq!(counts.iter().count(), 0)
}

#[test]
fn test_reverse_iteration_over_empty() {
    let counts = SampleAlleleCounts::default();
    assert_eq!(counts.iter().rev().count(), 0)
}

#[test]
fn test_count() {
    let mut counts = SampleAlleleCounts::default();
    counts.add_site_from_counts(AlleleCounts::try_new(&[1, 2, 3], 6).unwrap());
    assert_eq!(counts.iter().count(), 1);
}

#[cfg(test)]
fn make_nonempty_counts() -> SampleAlleleCounts {
    let mut counts = SampleAlleleCounts::default();
    counts.add_site_from_counts(AlleleCounts::try_new(&[1, 2, 3], 6).unwrap());
    counts.add_site_from_counts(AlleleCounts::try_new(&[1, 1, 1], 3).unwrap());
    counts.add_site_from_counts(AlleleCounts::try_new(&[1, 5, 1], 7).unwrap());
    counts.add_site_from_counts(AlleleCounts::try_new(&[1, 6, 2], 9).unwrap());

    counts
}

#[test]
fn test_iter_count() {
    let counts = make_nonempty_counts();
    assert_eq!(counts.iter().count(), counts.len());
    assert_eq!(counts.iter().filter(|c| c.counts()[1] == 5).count(), 1);

    let mut iter = counts.iter();
    let _ = iter.next().unwrap();
    assert_eq!(iter.count(), 3);
}

#[test]
fn test_nth() {
    let counts = make_nonempty_counts();
    let mut iter = counts.iter();
    assert_eq!(iter.nth(2), counts.get(2));
    let mut iter = counts.iter();
    let _ = iter.next().unwrap();
    assert_eq!(iter.nth(1), counts.get(2));
}

#[test]
fn test_nth_back() {
    let counts = make_nonempty_counts();
    let mut iter = counts.iter();
    assert_eq!(iter.nth_back(0), counts.get(3));
    assert_eq!(iter.nth_back(2), counts.get(0));
}

#[test]
fn test_exhaust_back() {
    // make sure we don't panic on decrementing 0usize
    let counts = make_nonempty_counts();
    let mut iter = counts.iter();
    _ = iter.next_back();
    _ = iter.nth_back(1);
    assert!(iter.next_back().is_some());
    assert!(iter.next_back().is_none());
}

#[test]
fn test_single_site_getters() {
    let mut counts = SampleAlleleCounts::default();
    counts.add_site_from_counts(AlleleCounts::try_new(&[9, 8, 7], 35).unwrap());
    let site = counts.get(0).unwrap();
    assert_eq!(site.counts(), &[9, 8, 7]);
    assert_eq!(site.total_alleles(), 35);
}
