use crate::iter::MultiSiteCountsIter;
#[cfg(feature = "tskit")]
use crate::{from_tree_sequence, from_tskit::FromTreeSequenceOptions};
use crate::{AlleleID, PopgenError, PopgenResult};
use std::cmp::max;

pub type Count = i64;

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
    /// Convenience wrapper which repeatedly invokes [`Self::add_site`].
    /// # Errors
    /// The error conditions from [`Self::add_site`] apply here.
    pub fn try_from_tabular<Sites, Samples>(sites: Sites) -> PopgenResult<Self>
    where
        Sites: IntoIterator<Item = Samples>,
        Samples: IntoIterator<Item = Option<AlleleID>>,
    {
        let mut ret = Self::default();

        for site in sites {
            ret.add_site(site)?;
        }

        Ok(ret)
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
    pub fn try_from_tree_sequence<N>(
        ts: &tskit::TreeSequence,
        samples: N,
        options: Option<FromTreeSequenceOptions>,
    ) -> Result<Self, PopgenError>
    where
        N: Iterator<Item = tskit::NodeId>,
    {
        Self::try_from_tree_sequence_site_iter(ts, samples, ts.site_iter(), options)
    }

    #[cfg(feature = "tskit")]
    pub fn try_from_tree_sequence_site_iter<'ts, N, S>(
        ts: &'ts tskit::TreeSequence,
        samples: N,
        sites: S,
        options: Option<FromTreeSequenceOptions>,
    ) -> Result<Self, PopgenError>
    where
        N: Iterator<Item = tskit::NodeId>,
        S: Iterator<Item = tskit::SiteRef<'ts>>,
    {
        from_tree_sequence::try_from_tree_sequence_with_site_iter(ts, samples, sites, options)
    }

    /// Add a site from an iterator of potentially missing allele IDs.
    /// # Errors
    /// - If `samples` is empty.
    /// - If `samples` contains no present data (i.e. only ever yields `None`).
    pub fn add_site<Samples>(&mut self, samples: Samples) -> PopgenResult<()>
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

        // we should not get NegativeCount or TotalAllelesDeficient here (could check that),
        // but we certainly could get other error variants
        self.add_site_from_counts(counts_this_site, total_alleles)
    }

    /// Add a site with counts of present alleles as described by `counts` and `total_alleles`
    /// alleles in total, including missing data.
    ///
    /// # Errors
    /// - Passed site counts must be valid.
    ///   See [`SiteCounts::try_new`].
    ///
    /// If any error occurs, the underlying struct has not been modified.
    // NOTE: any pub function that calls this one (probably)
    // should be fallible
    pub fn add_site_from_counts(
        &mut self,
        counts: impl AsRef<[Count]>,
        total_alleles: i32,
    ) -> PopgenResult<()> {
        let counts = counts.as_ref();
        if counts.is_empty() || total_alleles == 0 {
            return Err(PopgenError::EmptySiteCounts);
        }

        // check the conditions
        let _ = SiteCounts::try_new(counts, total_alleles)?;

        self.add_site_from_counts_unchecked(counts, total_alleles);
        Ok(())
    }

    fn add_site_from_counts_unchecked(&mut self, counts: &[Count], total_alleles: i32) {
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

/// A borrowed collection of allele counts and the total number of alleles (to describe, by implication, number of missing alleles).
///
/// This type is returned when requesting views into [`MultiSiteCounts`] and [`MultiPopulationCounts`].
/// It can also be built from user-provided data via [`Self::try_new`].
#[derive(Eq, PartialEq, Debug, Clone)]
pub struct SiteCounts<'inner> {
    pub(crate) counts: &'inner [Count],
    pub(crate) total_alleles: i32,
}

impl<'inner> SiteCounts<'inner> {
    /// Build a new `Self`, viewing a slice of counts and total alleles provided by the user.
    ///
    /// # Errors
    /// - If any element in `counts` is negative.
    /// - If `total_alleles` is less than the sum of elements of `counts`.
    /// - If `counts` is empty.
    /// - If `total_alleles == 0`.
    pub fn try_new(counts: &'inner [Count], total_alleles: i32) -> Result<Self, PopgenError> {
        if counts.is_empty() || total_alleles == 0 {
            return Err(PopgenError::EmptySiteCounts);
        }

        let mut sum = 0;
        for c in counts {
            if c < &0 {
                return Err(PopgenError::NegativeCount(*c));
            }
            sum += c;
        }

        if sum as i32 > total_alleles {
            return Err(PopgenError::TotalAllelesDeficient);
        }

        Ok(Self {
            counts,
            total_alleles,
        })
    }

    pub fn counts(&self) -> &[Count] {
        self.counts
    }

    pub fn total_alleles(&self) -> i32 {
        self.total_alleles
    }
}

/// Counts of present allele variants and of all alleles, including missing ones.
/// The data layout is by site, then population.
///
/// It is guaranteed that counts for the same site in multiple populations are meaningfully related,
/// particularly, e.g., that the allele assigned ID 0 in one population has also been assigned ID 0 in another population.
/// Alleles which appear in one population but not the other will have a count of 0 in that other population.
#[derive(Debug, Default, Clone)]
pub struct MultiPopulationCounts {
    // positions: Vec<usize>
    // ragged array (site, population) -> some collection of counts
    counts: Vec<Count>,
    // shape (site, population) -> index into counts
    count_starts: Vec<usize>,
    // (site, population) -> number of alleles, present or missing, at this site
    total_alleles: Vec<i32>,
    num_populations: usize,
}

impl MultiPopulationCounts {
    /// Create a new [`Self`] containing `how_many` populations, but containing no data.
    pub fn of_empty_populations(how_many: usize) -> Self {
        Self {
            counts: vec![],
            count_starts: vec![],
            total_alleles: vec![],
            num_populations: how_many,
        }
    }

    /// Extend the populations contained in [`Self`], using the successive (allele counts, number of samples) pairs provided.
    /// The first pair will be used to form the counts at this new site in the first population.
    /// The second pair will form the counts at the same site in the second population, etc.
    ///
    /// # Errors
    /// This function will fail **without rollback guarantees** if the provided slices do not match in length.
    /// The sites must also be individually valid; see [`SiteCounts::try_new`].
    /// Failure due to an invalid site does not provide rollback guarantees.
    pub fn extend_populations_from_site<Counts>(
        &mut self,
        mut get_counts: impl FnMut(usize) -> (Counts, i32),
    ) -> PopgenResult<()>
    where
        Counts: AsRef<[Count]>,
    {
        let mut inferred_slice_length = None;
        self.count_starts.push(self.counts.len());

        for population_i in 0..self.num_populations {
            let (allele_counts, total_alleles) = get_counts(population_i);
            let counts = allele_counts.as_ref();
            match inferred_slice_length {
                None => {
                    inferred_slice_length = Some(counts.len());
                    self.counts
                        .reserve(self.num_populations * inferred_slice_length.unwrap_or_default());
                }
                Some(stored) if stored != counts.len() => {
                    return Err(PopgenError::MismatchedSliceLength);
                }
                Some(_) => {}
            }

            let _ = SiteCounts::try_new(counts, total_alleles)?;

            self.counts.extend(counts);
            self.total_alleles.push(total_alleles);
        }

        Ok(())
    }

    #[cfg(feature = "tskit")]
    pub fn try_from_tree_sequence<Outer, Inner>(
        ts: &tskit::TreeSequence,
        samples: Outer,
        options: Option<FromTreeSequenceOptions>,
    ) -> Result<Self, PopgenError>
    where
        Outer: Iterator<Item = Inner>,
        Inner: Iterator<Item = tskit::NodeId>,
    {
        Self::try_from_tree_sequence_site_iter(ts, samples, ts.site_iter(), options)
    }

    #[cfg(feature = "tskit")]
    pub fn try_from_tree_sequence_site_iter<'ts, Outer, Inner, S>(
        ts: &'ts tskit::TreeSequence,
        samples: Outer,
        sites: S,
        options: Option<FromTreeSequenceOptions>,
    ) -> Result<Self, PopgenError>
    where
        Outer: Iterator<Item = Inner>,
        Inner: Iterator<Item = tskit::NodeId>,
        S: Iterator<Item = tskit::SiteRef<'ts>>,
    {
        from_tree_sequence::try_from_tree_sequence_multi_with_site_iter(ts, samples, sites, options)
    }

    /// Return the number of populations contained in [`Self`].
    pub fn num_populations(&self) -> usize {
        self.num_populations
    }

    /// Return the number of sites contained in `Self`.
    pub fn num_sites(&self) -> usize {
        self.total_alleles
            .len()
            .checked_div(self.num_populations)
            .unwrap_or(0)
    }

    /// Attempt to get a [`SiteCounts`] from `Self`.
    ///
    /// # Errors
    /// If any index is out of bounds.
    pub fn get(&self, site_num: usize, population_num: usize) -> Option<SiteCounts<'_>> {
        let counts_start = *self.count_starts.get(site_num)?;
        let counts_all_pops = match self.count_starts.get(site_num + 1) {
            None => &self.counts[counts_start..],
            Some(&next) => &self.counts[counts_start..next],
        };
        let counts_per_site = counts_all_pops.len().checked_div(self.num_populations())?;
        let counts_this_pop = &counts_all_pops
            [counts_per_site * population_num..counts_per_site * population_num + counts_per_site];

        let total_alleles = *self
            .total_alleles
            .get(site_num * population_num + population_num)?;

        Some(SiteCounts {
            counts: counts_this_pop,
            total_alleles,
        })
    }

    /// Convenience method to produce an iterator over [`SiteCounts`] in one population.
    /// Empty if `population_num` is out of bounds.
    pub fn iter_sites_in(&self, population_num: usize) -> impl Iterator<Item = SiteCounts<'_>> {
        (0..self.num_sites()).flat_map(move |site_n| self.get(site_n, population_num))
    }
}
