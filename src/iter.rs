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
}

impl DoubleEndedIterator for MultiSiteCountsIter<'_> {
    fn next_back(&mut self) -> Option<Self::Item> {
        let ret = self.inner.get(self.next_site_ind.1)?;

        self.next_site_ind.1 -= 1;
        Some(ret)
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
fn test_single_site_getters() {
    let mut counts = crate::MultiSiteCounts::default();
    counts.add_site_from_counts([9, 8, 7], 35);
    let site = counts.get(0).unwrap();
    assert_eq!(site.counts(), &[9, 8, 7]);
    assert_eq!(site.total_alleles(), 35);
}
