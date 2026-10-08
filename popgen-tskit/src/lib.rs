//! API to generate [`popgen::SampleAlleleCounts`] from [`tskit::TreeSequence`]

mod details;

use popgen::SampleAlleleCounts;

/// Options affecting the behavior of
/// functions processing tree sequence input.
#[derive(Debug, Default)]
pub struct FromTreeSequenceOptions {}

/// Error type related to tree sequence input
#[non_exhaustive]
#[derive(Debug)]
pub enum FromTreeSequenceError {
    /// Holds [`tskit::TskitError`]
    Tskit(::tskit::TskitError),
    /// A [`tskit::NodeId`] that is out of range with
    /// respect to the node table of a given tree sequence.
    NodeIdOutOfRange {
        /// The specific value that is out of range
        which: tskit::NodeId,
    },
    /// A site requires an ancestral state but none was present
    SiteMissingAncestralState,
    /// A mutation requires a derived state but none was present
    MutationMissingDerivedState,
    /// Position values were not sorted in increasing order
    UnsortedPositions,
    /// List of genomic windows is empty
    EmptyWindows,
    /// Returned if a window constains invalid positions,
    /// is not a proper interval,
    /// or if it overlaps with another window
    InvalidWindow((tskit::Position, tskit::Position)),
    /// Contains [`popgen::PopgenError`]
    Popgen(popgen::PopgenError),
}

/// Obtain site counts from a [`tskit::TreeSequence`].
/// All sites will be placed in one sample set.
/// If that is not desired, use [`multi_sample_allele_counts`] and related functions.
///
/// # Parameters
///
/// * `ts`: [`tskit::TreeSequence`]
/// * `samples`: iterator over [`tskit::NodeId`]
/// * `options`: modify the behavior using  [`FromTreeSequenceOptions`]
///
/// # Errors
///
/// Any errors from [`tskit`] will be propagated.
/// Missing ancestral/derived states will result in errors.
pub fn single_sample_allele_counts<N>(
    ts: &tskit::TreeSequence,
    samples: N,
    options: Option<FromTreeSequenceOptions>,
) -> Result<popgen::SampleAlleleCounts, FromTreeSequenceError>
where
    N: Iterator<Item = tskit::NodeId>,
{
    details::try_from_tree_sequence(ts, samples, options)
}

/// Obtain site counts from a [`tskit::TreeSequence`] using an iterator
/// over sites to include/exclude sites as needed.
///
/// See [`single_sample_allele_counts`] for details.
///
/// # Parameters
///
/// The additional parameter is:
///
/// `sites`: Iterator over [`tskit::SiteRef`]
pub fn single_sample_allele_counts_with_site_iter<'ts, N, S>(
    ts: &'ts tskit::TreeSequence,
    samples: N,
    sites: S,
    options: Option<FromTreeSequenceOptions>,
) -> Result<popgen::SampleAlleleCounts, FromTreeSequenceError>
where
    N: Iterator<Item = tskit::NodeId>,
    S: Iterator<Item = tskit::SiteRef<'ts>>,
{
    details::try_from_tree_sequence_with_site_iter(ts, samples, sites, options)
}

/// Obtain site counts from a [`tskit::TreeSequence`] using an iterator
/// over windows.
///
/// See [`single_sample_allele_counts`] for details.
/// API differences from that function are detailed below.
///
/// # Parameters
///
/// The additional parameter is:
///
/// `windows`: Iterator over a tuple of two elements convertible to
///            [`tskit::Position`] via [`From`].
///
/// # Returns
///
/// Unlike [`single_sample_allele_counts`], this function
/// returns a vector of [`popgen::SampleAlleleCounts`] with one element
/// per window.
///
/// # Notes
///
/// The tuples in `windows` represent half-open intervals, `[left, right)`.
pub fn single_sample_allele_counts_from_windows<N, W, P>(
    ts: &tskit::TreeSequence,
    samples: N,
    windows: W,
    options: Option<FromTreeSequenceOptions>,
) -> Result<Vec<popgen::SampleAlleleCounts>, FromTreeSequenceError>
where
    N: Iterator<Item = tskit::NodeId>,
    W: Iterator<Item = (P, P)>,
    P: Into<tskit::Position>,
{
    details::try_from_tree_sequence_windows(ts, samples, windows, options)
}

/// Construct count data from a tree sequence with respect to multiple
/// sample sets.
///
/// # Paramters
///
/// `ts`: [`tskit::TreeSequence`]
/// `samples`: Iterator over iterators of [`tskit::NodeId`]
/// `options`: [`FromTreeSequenceOptions`]
///
/// # Errors
///
/// Any errors from [`tskit`] will be propagated.
/// Missing ancestral/derived states will result in errors.
pub fn multi_sample_allele_counts<Outer, Inner>(
    ts: &tskit::TreeSequence,
    samples: Outer,
    options: Option<FromTreeSequenceOptions>,
) -> Result<popgen::SampleAlleleCounts, FromTreeSequenceError>
where
    Outer: Iterator<Item = Inner>,
    Inner: Iterator<Item = tskit::NodeId>,
{
    details::try_multi_sample_set_from_tree_sequence(ts, samples, options)
}

/// Construct count data from a tree sequence with respect to multiple
/// sample sets and site position ranges.
///
/// # Paramters
///
/// `ts`: [`tskit::TreeSequence`]
/// `samples`: Iterator over iterators of [`tskit::NodeId`]
/// `sites`: Iterator over [`tskit::SiteRef`]
/// `options`: [`FromTreeSequenceOptions`]
///
/// # Notes
///
/// The `sites` iterator should be obtained via [`tskit::TreeSequence::site_iter`]
/// and can be filtered using the [`Iterator`] API.
///
/// # Errors
///
/// Any errors from [`tskit`] will be propagated.
/// Missing ancestral/derived states will result in errors.
pub fn multi_sample_allele_counts_with_site_iter<'ts, Outer, Inner, S>(
    ts: &'ts tskit::TreeSequence,
    samples: Outer,
    sites: S,
    options: Option<FromTreeSequenceOptions>,
) -> Result<crate::SampleAlleleCounts, FromTreeSequenceError>
where
    Outer: Iterator<Item = Inner>,
    Inner: Iterator<Item = tskit::NodeId>,
    S: Iterator<Item = tskit::SiteRef<'ts>>,
{
    details::try_from_tree_sequence_multi_with_site_iter(ts, samples, sites, options)
}

impl std::fmt::Display for FromTreeSequenceError {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            FromTreeSequenceError::Tskit(e) => write!(f, "tskit error: {e}"),
            FromTreeSequenceError::NodeIdOutOfRange { which } => {
                write!(f, "node id {which} out of range")
            }
            FromTreeSequenceError::SiteMissingAncestralState => {
                write!(f, "site is missing ancestral state")
            }
            FromTreeSequenceError::MutationMissingDerivedState => {
                write!(f, "mutation is missing derived state")
            }
            FromTreeSequenceError::UnsortedPositions => {
                write!(f, "positions are not in increasing order")
            }
            FromTreeSequenceError::EmptyWindows => {
                write!(f, "empty windows")
            }
            FromTreeSequenceError::InvalidWindow(w) => {
                write!(f, "invalid window: {w:?}")
            }
            FromTreeSequenceError::Popgen(e) => {
                write!(f, "{e:?}")
            }
        }
    }
}

impl From<popgen::PopgenError> for FromTreeSequenceError {
    fn from(value: popgen::PopgenError) -> Self {
        Self::Popgen(value)
    }
}
