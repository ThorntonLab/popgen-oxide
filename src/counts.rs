use crate::iter::{SampleAlleleCountsSampleSetIter, SampleAlleleCountsSiteIter};
use crate::traits::TryReduce;
#[cfg(feature = "tskit")]
use crate::{from_tree_sequence, from_tskit::FromTreeSequenceOptions};
use crate::{AlleleID, PopgenError, PopgenResult};
use std::cmp::max;

/// A number of alleles.
/// Expected to be non-negative.
pub type Count = i64;

/// Counts of present allele variants and of all alleles, including missing ones.
/// The data layout is by site, then sample set.
///
/// When sample set numbers are used as references into this type, they are 0-based.
///
/// It is guaranteed that counts for the same site in multiple sample set are meaningfully related,
/// particularly, e.g., that the allele assigned ID 0 in one sample set has also been assigned ID 0 in another sample set.
/// Alleles which appear in one sample set but not the other will have a count of 0 in that other sample set, i.e. the counts will be padded to enforce the correspondence across sample sets.
///
/// The [`Default`] implementation currently creates a type with 0 sample sets.
/// This may not be what you want; consider calling [`Self::of_empty_sample_sets`].
#[derive(Debug, Default, Clone)]
pub struct SampleAlleleCounts {
    // probably don't need to track this
    // positions: Vec<i64>,
    counts: Vec<Count>,
    // start indices into counts at which the counts start for a specific site
    // counts and count_starts together produce a ragged 2d array
    count_starts: Vec<usize>,
    // (site, sample set) -> number of alleles, present or missing, at this site
    total_alleles: Vec<Count>,
    num_sample_sets: usize,
}

impl SampleAlleleCounts {
    /// Convenience wrapper which repeatedly invokes [`Self::add_site`].
    /// This will produce a `Self` with 1 sample set.
    ///
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

        // if no sites were added, this has not been set
        ret.num_sample_sets = 1;

        Ok(ret)
    }

    /// Obtain site counts from a [`tskit::TreeSequence`].
    /// All sites will be placed in one sample set; if that is not desired, use [`Self::try_multi_sample_set_from_tree_sequence`] and related functions.
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

    /// [`Self::try_from_tree_sequence`], but specifying a selection of sites using `sites`.
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

    /// [`Self::try_from_tree_sequence`], but specifying a selection of genomic windows using `windows`.
    ///
    /// Each window will be placed in a new `Self`, so this function returns a [`Vec`].
    #[cfg(feature = "tskit")]
    pub fn try_from_tree_sequence_windows<N, W, P>(
        ts: &tskit::TreeSequence,
        samples: N,
        windows: W,
        options: Option<FromTreeSequenceOptions>,
    ) -> Result<Vec<Self>, PopgenError>
    where
        N: Iterator<Item = tskit::NodeId>,
        W: Iterator<Item = (P, P)>,
        P: Into<tskit::Position>,
    {
        crate::from_tree_sequence::try_from_tree_sequence_windows(ts, samples, windows, options)
    }

    /// Add a site from an iterator of potentially missing allele IDs.
    ///
    /// It is assumed that `self` contains 0 or 1 sample sets.
    ///
    /// # Errors
    /// - If `self` contains more than 1 sample set.
    /// - If `samples` is empty.
    /// - If `samples` contains no present data (i.e. only ever yields `None`).
    pub fn add_site<Samples>(&mut self, samples: Samples) -> PopgenResult<()>
    where
        Samples: IntoIterator<Item = Option<AlleleID>>,
    {
        // TODO: multi sample set impl
        match self.num_sample_sets {
            0 => {
                self.num_sample_sets = 1;
            }
            1 => {}
            _more => {
                return Err(PopgenError::LibraryError(String::from(
                    "cannot add_site with more than one sample set",
                )));
            }
        }

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
        self.extend_sample_sets_from_site(|_| (&counts_this_site, total_alleles))?;
        Ok(())
    }

    /// Create a new [`Self`] containing `how_many` sample sets, but containing no data.
    pub fn of_empty_sample_sets(how_many: usize) -> Self {
        Self {
            counts: vec![],
            count_starts: vec![],
            total_alleles: vec![],
            num_sample_sets: how_many,
        }
    }

    /// Extend the sample sets contained in [`Self`], using the successive (allele counts, number of samples) pairs provided.
    /// The first pair will be used to form the counts at this new site in the first sample set.
    /// The second pair will form the counts at the same site in the second sample set, etc.
    ///
    /// Because of the invariant of this type, the same position in each counts slice must correspond to the same allele.
    /// Padding with zeroes may be needed to achieve this.
    ///
    /// # Errors
    /// This function will fail **without rollback guarantees** if the provided counts slices do not match in length.
    /// The sites must also be individually valid; see [`AlleleCounts::try_new`].
    /// Failure does not provide rollback guarantees.
    pub fn extend_sample_sets_from_site<Counts>(
        &mut self,
        mut get_counts: impl FnMut(usize) -> (Counts, Count),
    ) -> PopgenResult<()>
    where
        Counts: AsRef<[Count]>,
    {
        let mut inferred_slice_length = None;
        self.count_starts.push(self.counts.len());

        for sample_set_i in 0..self.num_sample_sets {
            let (allele_counts, total_alleles) = get_counts(sample_set_i);
            let counts = allele_counts.as_ref();
            match inferred_slice_length {
                None => {
                    inferred_slice_length = Some(counts.len());
                    self.counts
                        .reserve(self.num_sample_sets * inferred_slice_length.unwrap_or_default());
                }
                Some(stored) if stored != counts.len() => {
                    return Err(PopgenError::MismatchedSliceLength);
                }
                Some(_) => {}
            }

            let _ = AlleleCounts::try_new(counts, total_alleles)?;

            self.counts.extend(counts);
            self.total_alleles.push(total_alleles);
        }

        Ok(())
    }

    /// [`Self::try_from_tree_sequence`], with the ability to specify
    #[cfg(feature = "tskit")]
    pub fn try_multi_sample_set_from_tree_sequence<Outer, Inner>(
        ts: &tskit::TreeSequence,
        samples: Outer,
        options: Option<FromTreeSequenceOptions>,
    ) -> Result<Self, PopgenError>
    where
        Outer: Iterator<Item = Inner>,
        Inner: Iterator<Item = tskit::NodeId>,
    {
        Self::try_multi_sample_set_from_tree_sequence_site_iter(
            ts,
            samples,
            ts.site_iter(),
            options,
        )
    }

    #[cfg(feature = "tskit")]
    pub fn try_multi_sample_set_from_tree_sequence_site_iter<'ts, Outer, Inner, S>(
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

    /// Return the number of sample sets contained in [`Self`].
    pub fn num_sample_sets(&self) -> usize {
        self.num_sample_sets
    }

    /// `true` if and only if there are no sites in this [`Self`].
    pub fn is_empty(&self) -> bool {
        self.total_alleles.is_empty()
    }

    /// Return the number of sites contained in `Self`.
    pub fn num_sites(&self) -> usize {
        self.total_alleles
            .len()
            .checked_div(self.num_sample_sets)
            .unwrap_or(0)
    }

    /// Attempt to get a [`AlleleCounts`] from `Self` with respect to a given site and sample set.
    ///
    /// # Errors
    /// If any index is out of bounds.
    pub fn get_site(&self, site_num: usize, sample_set_num: usize) -> Option<AlleleCounts<'_>> {
        let counts_start = *self.count_starts.get(site_num)?;
        let counts_all_pops = match self.count_starts.get(site_num + 1) {
            None => &self.counts[counts_start..],
            Some(&next) => &self.counts[counts_start..next],
        };
        let counts_per_site = counts_all_pops.len().checked_div(self.num_sample_sets())?;
        let counts_this_pop = &counts_all_pops
            [counts_per_site * sample_set_num..counts_per_site * sample_set_num + counts_per_site];

        let total_alleles = *self
            .total_alleles
            .get(site_num * self.num_sample_sets() + sample_set_num)?;

        Some(AlleleCounts {
            counts: counts_this_pop,
            total_alleles,
        })
    }

    /// Convenience method equivalent to calling [`Self::iter_sample_set`] for each sample set in order.
    pub fn iter_sample_sets(
        &'_ self,
    ) -> impl DoubleEndedIterator<Item = SingleSampleAlleleCounts<'_>>
           + DoubleEndedIterator
           + ExactSizeIterator {
        SampleAlleleCountsSampleSetIter {
            inner: self,
            next_sample_set_ind: (0, self.num_sample_sets().saturating_sub(1)),
        }
    }

    /// Create the view [`SingleSampleAlleleCounts`] for the given sample set number.
    ///
    /// `None` if `sample_set_number` is out of bounds.
    pub fn sample_set(&'_ self, sample_set_number: usize) -> Option<SingleSampleAlleleCounts<'_>> {
        (0..self.num_sample_sets())
            .contains(&sample_set_number)
            .then_some(SingleSampleAlleleCounts {
                inner: self,
                sample_set_number,
            })
    }

    /// Shortcut via [`Self::sample_set`] to [`SingleSampleAlleleCounts::into_iter`].
    pub fn iter_sample_set(
        &'_ self,
        sample_set_number: usize,
    ) -> Option<SampleAlleleCountsSiteIter<'_>> {
        Some(self.sample_set(sample_set_number)?.into_iter())
    }
}

impl TryReduce for SampleAlleleCounts {
    type Error = PopgenError;

    /// Attempt to concatenate `self` and `other`, assuming that the sample sets correspond, with the semantics that the sites from `self` will be followed by the sites from `other`.
    ///
    /// Error if the number of sample sets differs.
    fn try_reduce(self, other: Self) -> Result<Self, Self::Error>
    where
        Self: Sized,
    {
        if self.num_sample_sets != other.num_sample_sets {
            return Err(PopgenError::MismatchedSampleSetCount(
                self.num_sample_sets(),
                other.num_sample_sets(),
            ));
        }

        let counts_len_left = self.counts.len();
        let mut counts = self.counts;
        counts.extend(other.counts);

        let mut count_starts = self.count_starts;
        count_starts.extend(
            other
                .count_starts
                .into_iter()
                .map(|cs| counts_len_left + cs),
        );

        let mut total_alleles = self.total_alleles;
        total_alleles.extend(other.total_alleles);

        Ok(Self {
            counts,
            count_starts,
            total_alleles,
            num_sample_sets: self.num_sample_sets,
        })
    }
}

/// A view into [`SampleAlleleCounts`], equivalent to a reference to that type and a sample set number.
///
/// This struct may be consumed via [`IntoIterator`] to iterate over sites within this sample set.
#[derive(Debug, Clone)]
pub struct SingleSampleAlleleCounts<'i> {
    inner: &'i SampleAlleleCounts,
    sample_set_number: usize,
}

impl<'i> IntoIterator for SingleSampleAlleleCounts<'i> {
    type Item = AlleleCounts<'i>;
    type IntoIter = SampleAlleleCountsSiteIter<'i>;

    fn into_iter(self) -> Self::IntoIter {
        SampleAlleleCountsSiteIter {
            inner: self.inner,
            sample_set_number: self.sample_set_number,
            next_site_ind: (0, self.inner.num_sites().saturating_sub(1)),
        }
    }
}

impl<'i> SingleSampleAlleleCounts<'i> {
    /// Get the counts for this `site_number` within this sample set, erroring if out of range
    pub fn site(&'i self, site_number: usize) -> Option<AlleleCounts<'i>> {
        self.inner.get_site(site_number, self.sample_set_number)
    }

    pub fn inner(&self) -> &'i SampleAlleleCounts {
        self.inner
    }

    pub fn sample_set_number(&self) -> usize {
        self.sample_set_number
    }
}

/// A borrowed collection of allele counts and the total number of alleles (to describe, by implication, number of missing alleles).
///
/// This type is returned when requesting views into [`SampleAlleleCounts`].
/// It can also be built from user-provided data via [`Self::try_new`].
#[derive(Eq, PartialEq, Debug, Clone)]
pub struct AlleleCounts<'inner> {
    counts: &'inner [Count],
    total_alleles: i64,
}

impl<'inner> AlleleCounts<'inner> {
    /// Build a new `Self`, viewing a slice of counts and total alleles provided by the user.
    ///
    /// # Errors
    /// - If any element in `counts` is negative.
    /// - If `total_alleles` is less than the sum of elements of `counts`.
    /// - If `counts` is empty.
    /// - If `total_alleles == 0`.
    pub fn try_new(counts: &'inner [Count], total_alleles: i64) -> Result<Self, PopgenError> {
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

        if sum > total_alleles {
            return Err(PopgenError::TotalAllelesDeficient);
        }

        Ok(Self {
            counts,
            total_alleles,
        })
    }

    #[inline]
    pub fn counts(&self) -> &[Count] {
        self.counts
    }

    #[inline]
    pub fn total_alleles(&self) -> i64 {
        self.total_alleles
    }
}
