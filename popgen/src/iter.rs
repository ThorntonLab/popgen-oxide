use crate::{AlleleCounts, Count};

pub struct AlleleCountsSiteIter<'inner> {
    pub(crate) inner: &'inner AlleleCounts,
    // index of next for forward iter, index of next for reverse iter
    pub(crate) next_site_ind: (usize, usize),
}

pub struct IteratedSite<'inner> {
    counts: &'inner [Count],
    alleles_missing: i32,
}

impl<'inner> Iterator for AlleleCountsSiteIter<'inner> {
    type Item = IteratedSite<'inner>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.next_site_ind.0 > self.next_site_ind.1 {
            return None;
        }

        self.next_site_ind.0 += 1;
        Some(IteratedSite {
            counts: self.inner.counts_at(self.next_site_ind.0)
                .expect("forward iterator index out of range"),
            alleles_missing: self.inner.alleles_missing[self.next_site_ind.0],
        })
    }
}

impl DoubleEndedIterator for AlleleCountsSiteIter<'_> {
    fn next_back(&mut self) -> Option<Self::Item> {
        if self.next_site_ind.1 < self.next_site_ind.0 {
            return None;
        }

        self.next_site_ind.1 -= 1;
        Some(IteratedSite {
            counts: self.inner.counts_at(self.next_site_ind.1)
                .expect("reverse iterator index out of range"),
            alleles_missing: self.inner.alleles_missing[self.next_site_ind.1],
        })
    }
}