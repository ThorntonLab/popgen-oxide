use crate::iter::{MultiSiteCountsIter, SiteCounts};
use crate::{AlleleID, Count};
use std::cmp::max;
use std::ops::Index;

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

    /// Obtain site counts from a [`tskit::TreeSequence`].
    ///
    /// # Parameters
    ///
    /// * `ts`: [`tskit::TreeSequence`]
    /// * `options`: modify the behavior using  [`FromTreeSequenceOptions`]
    ///
    /// # Errors
    ///
    /// Any errors from [`tskit`] will be propagated.
    ///
    /// # Panics
    ///
    /// Sites with empty ancestral states and mutations with empty
    /// derived states are currently rejected as a hard error resulting
    /// in a panic.
    #[cfg(feature = "tskit")]
    pub fn try_from_tree_sequence(
        ts: &tskit::TreeSequence,
        options: Option<crate::FromTreeSequenceOptions>,
    ) -> Result<Self, crate::PopgenError> {
        crate::from_tree_sequence::try_from_tree_sequence(ts, options)
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

    /// The number of sites added to this [`Self`] so far.
    pub fn len(&self) -> usize {
        self.total_alleles.len()
    }

    /// `true` if and only if there are no sites in this [`Self`]; equivalent to `self.len() == 0`.
    pub fn is_empty(&self) -> bool {
        self.total_alleles.is_empty()
    }

    pub fn iter(&self) -> MultiSiteCountsIter<'_> {
        MultiSiteCountsIter {
            inner: self,
            next_site_ind: (
                0,
                if !self.count_starts.is_empty() {
                    self.count_starts.len() - 1
                } else {
                    0
                },
            ),
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
    pub fn get(&self, site: usize) -> Option<SiteCounts<'_>> {
        Some(SiteCounts {
            counts: self.counts_slice_at(site)?,
            total_alleles: self.total_alleles[site],
        })
    }
}

/// A collection of [`MultiSiteCounts`], with an added invariant.
///
/// A naive collection of such counts does not necessarily form a meaningful collection of populations with comparable counts.
/// This is because each implicitly encodes the mapping from actual allele to a position in the array describing the count of each allele.
/// This mapping may not be the same between independent [`MultiSiteCounts`] structs.
///
/// This struct remedies this by building these [`MultiSiteCounts`] such that they all respect the same mapping.
#[derive(Debug, Default, Clone)]
pub struct MultiPopulationCounts {
    pub(crate) populations: Vec<MultiSiteCounts>,
}

impl Index<usize> for MultiPopulationCounts {
    type Output = MultiSiteCounts;

    fn index(&self, index: usize) -> &Self::Output {
        &self.populations[index]
    }
}
