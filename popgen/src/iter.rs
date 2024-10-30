use crate::{Count, MultiSiteCounts};

pub struct MultiSiteCountsIter<'inner> {
    pub(crate) inner: &'inner MultiSiteCounts,
    // index of next for forward iter, index of next for reverse iter
    pub(crate) next_site_ind: (usize, usize),
}

#[derive(PartialEq, Debug)]
pub struct SiteCounts<'inner> {
    pub(crate) counts: &'inner [Count],
    pub(crate) alleles_missing: i32,
}

impl<'inner> Iterator for MultiSiteCountsIter<'inner> {
    type Item = SiteCounts<'inner>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.next_site_ind.0 > self.next_site_ind.1 {
            return None;
        }

        let ret = SiteCounts {
            counts: self.inner.counts_at(self.next_site_ind.0)
                .expect("forward iterator index out of range"),
            alleles_missing: self.inner.alleles_missing[self.next_site_ind.0],
        };

        self.next_site_ind.0 += 1;
        Some(ret)
    }
}

impl DoubleEndedIterator for MultiSiteCountsIter<'_> {
    fn next_back(&mut self) -> Option<Self::Item> {
        if self.next_site_ind.1 < self.next_site_ind.0 {
            return None;
        }

        let ret = SiteCounts {
            counts: self.inner.counts_at(self.next_site_ind.1)
                .expect("reverse iterator index out of range"),
            alleles_missing: self.inner.alleles_missing[self.next_site_ind.1],
        };

        self.next_site_ind.1 -= 1;
        Some(ret)
    }
}