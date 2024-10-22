use crate::iter::AlleleCountsSiteIter;
use std::cmp::max;
use std::error::Error;
use std::fmt::{Debug, Display, Formatter};
use std::str::FromStr;

mod tests;
mod iter;
pub mod adapter;

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

type Count = i64;

#[derive(Debug, Default)]
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

        self.add_site_from_counts(counts_this_site);
    }

    pub fn add_site_from_counts(&mut self, counts: impl AsRef<[Count]>) {
        self.counts.extend_from_slice(&*counts);
        // count backwards in case counts_this_site.is_empty() or other strange case
        self.count_starts.push(self.counts.len() - counts.len());
    }

    pub fn iter(&self) -> AlleleCountsSiteIter {
        AlleleCountsSiteIter {
            inner: &self,
            next_site_ind: (0, self.count_starts.len() - 1),
        }
    }

    /// Get the allele counts at a specific site index.
    ///
    /// Returns [`None`] if the site index is invalid.
    pub fn counts_at(&self, site: usize) -> Option<&[Count]> {
        self.count_starts.get(site).map(|count_start| &self.counts[*count_start..
            *self.count_starts
                .get(site + 1)
                .unwrap_or(&(self.counts.len() - 1))])
    }

    fn heterozyosity_from_slice(counts: &[Count]) -> f64 {
        let num_pairs = {
            let count: i64 = counts.iter().sum();
            count * (count - 1)
        };

        // the number of pairs where the two samples are homozygous, summed over every genotype
        let num_homozygous_pairs: Count = counts.iter().map(|count| count * (count - 1)).sum();

        1f64 - (num_homozygous_pairs as f64 / num_pairs as f64)
    }

    /// The chance pi that two uniformly randomly chosen genotypes at this site are different.
    /// Expect a value in the 1e-3 range.
    /// The complement of this probability is the site homozygosity.
    ///
    ///
    /// Returns [`None`] if the site index is invalid.
    pub fn site_heterozygosity(&self, site: usize) -> Option<f64> {
        self.counts_at(site).map(Self::heterozyosity_from_slice)
    }

    /// Take the sum of [`Self::site_heterozygosity`] over all sites.
    /// This statistic is the expected number of differences between the genotypes of two uniformly chosen individuals, considering all sites.
    pub fn global_heterozygosity(&self) -> f64 {
        self.iter().map(Self::heterozyosity_from_slice).sum()
    }
}

