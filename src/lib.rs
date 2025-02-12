use crate::iter::{MultiSiteCountsIter, SiteCounts};
use std::cmp::max;
use std::fmt::Debug;
use std::str::FromStr;

pub mod adapter;
pub mod iter;
pub mod stats;
mod test;
pub(crate) mod util;

pub use noodles::vcf as noodles_vcf;

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

// TODO: errors?

type Count = i64;

#[derive(Debug, Default, Clone)]
pub struct MultiSiteCounts {
    // probably don't need to track this
    // positions: Vec<i64>,
    counts: Vec<Count>,
    // start indices into counts at which the counts start for a specific site
    // counts and count_starts together produce a ragged 2d array
    count_starts: Vec<usize>,
    total_alleles: Vec<i32>,
}

impl MultiSiteCounts {
    pub fn from_tabular<Sites, Samples>(sites: Sites) -> Self
    where
        Sites: IntoIterator<Item = Samples>,
        Samples: IntoIterator<Item = Option<AlleleID>>,
    {
        let mut ret = Self::default();

        for site in sites {
            ret.add_site(site);
        }

        ret
    }

    pub fn add_site<Samples>(&mut self, samples: Samples)
    where
        Samples: IntoIterator<Item = Option<AlleleID>>,
    {
        let mut total_alleles = 0;

        // in something like VCF we wouldn't even have data if there was no variation; 2 is a reasonable lower bound
        // we're allocating `usize`s; it's totally fine to do this
        let mut counts_this_site = Vec::with_capacity(2);
        for allele_id in samples {
            total_alleles += 1;
            let allele_id_under = match allele_id {
                None => continue,
                Some(id) => id.0,
            };

            counts_this_site.resize(max(allele_id_under + 1, counts_this_site.len()), 0);
            counts_this_site[allele_id_under] += 1;
        }

        self.add_site_from_counts(counts_this_site, total_alleles);
    }

    pub fn add_site_from_counts(&mut self, counts: impl AsRef<[Count]>, total_alleles: i32) {
        self.counts.extend_from_slice(counts.as_ref());
        // count backwards in case counts_this_site.is_empty() or other strange case
        self.count_starts
            .push(self.counts.len() - counts.as_ref().len());
        self.total_alleles.push(total_alleles);
    }

    pub fn iter(&self) -> MultiSiteCountsIter {
        MultiSiteCountsIter {
            inner: self,
            next_site_ind: (0, self.count_starts.len() - 1),
        }
    }

    fn counts_slice_at(&self, site: usize) -> Option<&[Count]> {
        // TODO: do we even need random access?

        self.count_starts.get(site).map(|count_start| {
            &self.counts[*count_start
                ..*self
                    .count_starts
                    .get(site + 1)
                    .unwrap_or(&self.counts.len())]
        })
    }

    /// Get the allele counts at a specific site index.
    ///
    /// Returns [`None`] if the site index is invalid.
    pub fn counts_at(&self, site: usize) -> Option<SiteCounts> {
        Some(SiteCounts {
            counts: self.counts_slice_at(site)?,
            total_alleles: self.total_alleles[site],
        })
    }
}
