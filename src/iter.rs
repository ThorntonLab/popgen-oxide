use crate::{AlleleCounts, SampleAlleleCounts};

pub struct SampleAlleleCountsPopulationIter<'inner> {
    pub(crate) inner: &'inner SampleAlleleCounts,
    // index of next for forward iter, index of next for reverse iter
    pub(crate) next_population_ind: (usize, usize),
}

impl<'inner> Iterator for SampleAlleleCountsPopulationIter<'inner> {
    type Item = SampleAlleleCountsSiteIter<'inner>;

    fn next(&mut self) -> Option<Self::Item> {
        let ret = self.inner.iter_sites_in(self.next_population_ind.0)?;

        self.next_population_ind.0 += 1;
        Some(ret)
    }

    fn count(self) -> usize {
        self.len()
    }

    fn last(self) -> Option<Self::Item> {
        let mut s = self;
        s.next_population_ind.0 = s.next_population_ind.1;
        s.next()
    }

    // recall that skip uses this internally
    fn nth(&mut self, n: usize) -> Option<Self::Item> {
        self.next_population_ind.0 += n;
        self.next()
    }
}

impl ExactSizeIterator for SampleAlleleCountsPopulationIter<'_> {
    fn len(&self) -> usize {
        if self.inner.num_populations() == 0 {
            0
        } else {
            self.next_population_ind.1 - self.next_population_ind.0 + 1
        }
    }
}

impl DoubleEndedIterator for SampleAlleleCountsPopulationIter<'_> {
    fn next_back(&mut self) -> Option<Self::Item> {
        let ret = self.inner.iter_sites_in(self.next_population_ind.1)?;

        // there is no way to have usize counts because a Vec can never exceed isize::MAX
        self.next_population_ind.1 = self.next_population_ind.1.wrapping_sub(1);
        Some(ret)
    }

    fn nth_back(&mut self, n: usize) -> Option<Self::Item> {
        self.next_population_ind.1 = self.next_population_ind.1.wrapping_sub(n);
        self.next_back()
    }
}

pub struct SampleAlleleCountsSiteIter<'inner> {
    pub(crate) inner: &'inner SampleAlleleCounts,
    pub(crate) population_number: usize,
    // index of next for forward iter, index of next for reverse iter
    pub(crate) next_site_ind: (usize, usize),
}

impl<'inner> Iterator for SampleAlleleCountsSiteIter<'inner> {
    type Item = AlleleCounts<'inner>;

    fn next(&mut self) -> Option<Self::Item> {
        let ret = self
            .inner
            .get_site(self.next_site_ind.0, self.population_number)?;

        self.next_site_ind.0 += 1;
        Some(ret)
    }

    fn count(self) -> usize {
        self.len()
    }

    fn last(self) -> Option<Self::Item> {
        let mut s = self;
        s.next_site_ind.0 = s.next_site_ind.1;
        s.next()
    }

    // recall that skip uses this internally
    fn nth(&mut self, n: usize) -> Option<Self::Item> {
        self.next_site_ind.0 += n;
        self.next()
    }
}

impl ExactSizeIterator for SampleAlleleCountsSiteIter<'_> {
    fn len(&self) -> usize {
        if self.inner.is_empty() {
            0
        } else {
            self.next_site_ind.1 - self.next_site_ind.0 + 1
        }
    }
}

impl DoubleEndedIterator for SampleAlleleCountsSiteIter<'_> {
    fn next_back(&mut self) -> Option<Self::Item> {
        let ret = self
            .inner
            .get_site(self.next_site_ind.1, self.population_number)?;

        // there is no way to have usize counts because a Vec can never exceed isize::MAX
        self.next_site_ind.1 = self.next_site_ind.1.wrapping_sub(1);
        Some(ret)
    }

    fn nth_back(&mut self, n: usize) -> Option<Self::Item> {
        self.next_site_ind.1 = self.next_site_ind.1.wrapping_sub(n);
        self.next_back()
    }
}

#[test]
fn test_iteration_over_empty() {
    let counts = SampleAlleleCounts::default();
    assert_eq!(counts.iter_populations().count(), 0)
}

#[test]
fn test_reverse_iteration_over_empty() {
    let counts = SampleAlleleCounts::default();
    assert_eq!(counts.iter_populations().rev().count(), 0)
}

#[test]
fn test_population_count() {
    let mut counts = SampleAlleleCounts::of_empty_populations(1);
    counts
        .extend_populations_from_site(|_| (&[1, 2, 3], 6))
        .unwrap();

    assert_eq!(counts.iter_populations().count(), 1);
}

#[cfg(test)]
fn make_nonempty_counts() -> SampleAlleleCounts {
    let mut counts = SampleAlleleCounts::of_empty_populations(1);
    counts
        .extend_populations_from_site(|_| (&[1, 2, 3], 6))
        .unwrap();
    counts
        .extend_populations_from_site(|_| (&[1, 1, 1], 3))
        .unwrap();
    counts
        .extend_populations_from_site(|_| (&[1, 5, 1], 7))
        .unwrap();
    counts
        .extend_populations_from_site(|_| (&[1, 6, 2], 9))
        .unwrap();

    counts
}

#[test]
fn test_site_count() {
    let counts = make_nonempty_counts();
    assert_eq!(counts.iter_sites_in(0).unwrap().count(), counts.num_sites());
    assert_eq!(
        counts
            .iter_sites_in(0)
            .unwrap()
            .filter(|c| c.counts()[1] == 5)
            .count(),
        1
    );

    let mut iter = counts.iter_sites_in(0).unwrap();
    let _ = iter.next().unwrap();
    assert_eq!(iter.count(), 3);
}

#[test]
fn test_site_nth() {
    let counts = make_nonempty_counts();
    let mut iter = counts.iter_sites_in(0).unwrap();
    assert_eq!(iter.nth(2), counts.get_site(2, 0));
    let mut iter = counts.iter_sites_in(0).unwrap();
    let _ = iter.next().unwrap();
    assert_eq!(iter.nth(1), counts.get_site(2, 0));
}

#[test]
fn test_site_nth_back() {
    let counts = make_nonempty_counts();
    let mut iter = counts.iter_sites_in(0).unwrap();
    assert_eq!(iter.nth_back(0), counts.get_site(3, 0));
    assert_eq!(iter.nth_back(2), counts.get_site(0, 0));
}

#[test]
fn test_site_exhaust_back() {
    // make sure we don't panic on decrementing 0usize
    let counts = make_nonempty_counts();
    let mut iter = counts.iter_sites_in(0).unwrap();
    _ = iter.next_back();
    _ = iter.nth_back(1);
    assert!(iter.next_back().is_some());
    dbg!(iter.next_site_ind);
    assert!(iter.next_back().is_none());
}

#[test]
fn test_single_site_getters() {
    let mut counts = SampleAlleleCounts::of_empty_populations(1);
    counts
        .extend_populations_from_site(|_| (&[9, 8, 7], 35))
        .unwrap();
    let site = counts.get_site(0, 0).unwrap();
    assert_eq!(site.counts(), &[9, 8, 7]);
    assert_eq!(site.total_alleles(), 35);
}
