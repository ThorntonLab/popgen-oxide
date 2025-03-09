use crate::{Count, MultiSiteCounts};

pub struct MultiSiteCountsIter<'inner> {
    pub(crate) inner: &'inner MultiSiteCounts,
    // index of next for forward iter, index of next for reverse iter
    pub(crate) next_site_ind: (usize, usize),
}

#[derive(PartialEq, Debug, Clone)]
pub struct SiteCounts<'inner> {
    pub(crate) counts: &'inner [Count],
    pub(crate) total_alleles: i32,
}

impl SiteCounts<'_> {
    pub fn counts(&self) -> &[Count] {
        self.counts
    }

    pub fn total_alleles(&self) -> i32 {
        self.total_alleles
    }
}

impl<'inner> Iterator for MultiSiteCountsIter<'inner> {
    type Item = SiteCounts<'inner>;

    fn next(&mut self) -> Option<Self::Item> {
        let ret = self.inner.get(self.next_site_ind.0)?;

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

impl DoubleEndedIterator for MultiSiteCountsIter<'_> {
    fn next_back(&mut self) -> Option<Self::Item> {
        let ret = self.inner.get(self.next_site_ind.1)?;

        // there is no way to have usize counts because a Vec can never exceed isize::MAX
        self.next_site_ind.1 = self.next_site_ind.1.wrapping_sub(1);
        Some(ret)
    }

    fn nth_back(&mut self, n: usize) -> Option<Self::Item> {
        self.next_site_ind.1 = self.next_site_ind.1.wrapping_sub(n);
        self.next_back()
    }
}

impl ExactSizeIterator for MultiSiteCountsIter<'_> {
    fn len(&self) -> usize {
        if self.inner.is_empty() {
            0
        } else {
            self.next_site_ind.1 - self.next_site_ind.0 + 1
        }
    }
}

#[test]
fn test_iteration_over_empty() {
    let counts = crate::MultiSiteCounts::default();
    assert_eq!(counts.iter().count(), 0)
}

#[test]
fn test_reverse_iteration_over_empty() {
    let counts = crate::MultiSiteCounts::default();
    assert_eq!(counts.iter().rev().count(), 0)
}

#[test]
fn test_count() {
    let mut counts = crate::MultiSiteCounts::default();
    counts.add_site_from_counts([1, 2, 3], 6);
    assert_eq!(counts.iter().count(), 1);
}

#[cfg(test)]
fn make_nonempty_counts() -> crate::MultiSiteCounts {
    let mut counts = crate::MultiSiteCounts::default();
    counts.add_site_from_counts([1, 2, 3], 6);
    counts.add_site_from_counts([1, 1, 1], 3);
    counts.add_site_from_counts([1, 5, 1], 7);
    counts.add_site_from_counts([1, 6, 2], 9);

    counts
}

#[test]
fn test_iter_count() {
    let counts = make_nonempty_counts();
    assert_eq!(counts.iter().count(), counts.len());
    assert_eq!(counts.iter().filter(|c| c.counts[1] == 5).count(), 1);

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
    let mut counts = crate::MultiSiteCounts::default();
    counts.add_site_from_counts([9, 8, 7], 35);
    let site = counts.get(0).unwrap();
    assert_eq!(site.counts(), &[9, 8, 7]);
    assert_eq!(site.total_alleles(), 35);
}
