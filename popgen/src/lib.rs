use itertools::Itertools;
use std::cmp::max;
use std::error::Error;
use std::fmt::{Debug, Display, Formatter};
use std::str::FromStr;

mod tests;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct AlleleID(usize);

impl FromStr for AlleleID {
    type Err = <usize as FromStr>::Err;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Ok(Self(usize::from_str(s)?))
    }
}

impl From<usize> for AlleleID {
    fn from(value: usize) -> Self {
        Self(value)
    }
}


#[derive(Debug)]
struct PopgenError {
    msg: String,
}

impl Display for PopgenError {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(f, "{:?}", self.msg)
    }
}

impl Error for PopgenError {}

type Count = i32;

#[derive(Debug)]
pub struct AlleleCounts {
    // probably don't need to track this
    // positions: Vec<i64>,
    counts: Vec<Count>,
    // start indices into counts at which the counts start for a specific site
    // counts and count_starts together produce a ragged 2d array
    count_starts: Vec<usize>,
}

impl AlleleCounts {
    pub fn from_tabular<Sites, Samples>(sites: Sites) -> Self
    where
        Sites: IntoIterator<Item = Samples>,
        Samples: IntoIterator<Item = Option<AlleleID>>,
    {
        let mut ret = Self {
            counts: vec![],
            count_starts: vec![],
        };

        for site in sites {
            ret.add_site(site);
        }

        ret
    }

    pub fn add_site<Samples>(&mut self, samples: Samples)
    where
        Samples: IntoIterator<Item = Option<AlleleID>>,
    {
        // in something like VCF we wouldn't even have data if there was no variation; 2 is a reasonable lower bound
        // we're allocating `usize`s; it's totally fine to do this
        let mut counts_this_site = Vec::with_capacity(2);
        for allele_id in samples {
            let allele_id_under = match allele_id {
                None => continue,  // TODO: anything else to do here?
                Some(id) => id.0,
            };

            counts_this_site.resize(max(allele_id_under + 1, counts_this_site.len()), 0);
            counts_this_site[allele_id_under] += 1;
        }
        self.counts.extend_from_slice(&*counts_this_site);
        // count backwards in case counts_this_site.is_empty() or other strange case
        self.count_starts.push(self.counts.len() - counts_this_site.len());
    }

    pub fn iter(&self) -> AlleleCountsSiteIter {
        AlleleCountsSiteIter {
            inner: &self,
            next_site_ind: (0, self.count_starts.len() - 1),
        }
    }
}

pub struct AlleleCountsSiteIter<'counts> {
    inner: &'counts AlleleCounts,
    // index of next for forward iter, index of next for reverse iter
    next_site_ind: (usize, usize),
}

impl Iterator for AlleleCountsSiteIter<'counts> {
    type Item<'counts> = &'counts [Count];

    fn next(&mut self) -> Option<Self::Item> {
        if self.next_site_ind.0 > self.next_site_ind.1 {
            return None;
        }

        match self.inner.count_starts.get(self.next_site_ind.0) {
            None => None,
            Some(index) => {
                let ret = self.inner.counts[index..
                    self.inner.count_starts
                        .get(self.next_site_ind.0 + 1)
                        .unwrap_or(&(self.inner.counts.len() - 1))];

                self.next_site_ind.0 += 1;
                ret
            }
        }
    }
}

impl DoubleEndedIterator for AlleleCountsSiteIter {
    fn next_back(&mut self) -> Option<Self::Item> {
        if self.next_site_ind.1 < self.next_site_ind.0 {
            return None;
        }

        match self.inner.count_starts.get(self.next_site_ind.1) {
            None => None,
            Some(index) => {
                let ret = self.inner.counts[index..
                    self.inner.count_starts
                        .get(self.next_site_ind.1 + 1)
                        .unwrap_or(&(self.inner.counts.len() - 1))];

                self.next_site_ind.1 -= 1;
                ret
            }
        }
    }
}
