use crate::{AlleleCounts, Count};

pub struct AlleleCountsSiteIter<'inner> {
    pub(crate) inner: &'inner AlleleCounts,
    // index of next for forward iter, index of next for reverse iter
    pub(crate) next_site_ind: (usize, usize),
}

impl<'inner> Iterator for AlleleCountsSiteIter<'inner> {
    type Item = &'inner [Count];

    fn next(&mut self) -> Option<Self::Item> {
        if self.next_site_ind.0 > self.next_site_ind.1 {
            return None;
        }

        match self.inner.count_starts.get(self.next_site_ind.0) {
            None => None,
            Some(index) => {
                let ret = self.inner.counts_at(*index)
                    .expect("forward iterator index out of range");

                self.next_site_ind.0 += 1;
                Some(ret)
            }
        }
    }
}

impl DoubleEndedIterator for AlleleCountsSiteIter<'_> {
    fn next_back(&mut self) -> Option<Self::Item> {
        if self.next_site_ind.1 < self.next_site_ind.0 {
            return None;
        }

        match self.inner.count_starts.get(self.next_site_ind.1) {
            None => None,
            Some(index) => {
                let ret = self.inner.counts_at(*index)
                    .expect("reverse iterator index out of range");

                self.next_site_ind.1 -= 1;
                Some(ret)
            }
        }
    }
}