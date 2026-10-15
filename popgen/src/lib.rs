#![cfg_attr(doc_cfg, feature(doc_cfg))]
//! Efficient types and interfaces for population genetics
//!
//! # Examples
//!
//! ## Collecting allele count data into [SampleAlleleCounts] for a single sample set.
//!
//! ```rust
//! let mut counts =
//! popgen::SampleAlleleCounts::of_empty_sample_sets(std::num::NonZeroUsize::new(1).unwrap());
//! let c = vec![10_i64, 2];
//! // NOTE: the number 14 is larger than the sum of elements in c,
//! // which means that there are 2 missing allele calls in the data.
//! counts.extend_sample_sets_from_site_iter(std::iter::once((&c, 14))).unwrap();
//! assert_eq!(counts.num_sites(), 1);
//! ```
//! ## Calculating diversity from [SampleAlleleCounts]
//!
//! ```rust
//! # let mut counts =
//! popgen::SampleAlleleCounts::of_empty_sample_sets(std::num::NonZeroUsize::new(1).unwrap());
//! # let c = vec![10_i64, 2];
//! # // NOTE: the number 14 is larger than the sum of elements in c,
//! # // which means that there are 2 missing allele calls in the data.
//! # counts.extend_sample_sets_from_site_iter(std::iter::once((&c, 14))).unwrap();
//! use popgen::stats::UnpolarisedSiteStat;
//! let d = popgen::stats::Diversity::try_from_iter_sites(counts.iter_sample_set(0).unwrap()).unwrap();
//! let manual = 1.0 - c.iter().map(|&i| (i as f64)*((i-1) as f64)).sum::<f64>()/(12.*11.);
//! assert_eq!(d, manual);
//! let d_f64 = f64::from(d);
//! assert_eq!(d_f64, manual);
//! // more rigorous method:
//! assert!((d_f64 - manual).abs() <= 1e-10);
//! ```
//!
//! ## Calculating diversity from an iterator over slices of count data.
//!
//! ```rust
//! use popgen::stats::UnpolarisedSiteStat;
//! use popgen::AlleleCounts;
//! let c = vec![vec![10_i64, 2], vec![8, 6], vec![1, 6, 3]];
//! let iter = c.iter().map(|v| AlleleCounts::new(v.as_slice(), 14).unwrap());
//! let d = popgen::stats::Diversity::try_from_iter_sites(iter).unwrap();
//! ```

use std::fmt::Debug;

mod counts;
pub mod iter;
pub mod stats;
#[allow(missing_docs)]
pub mod traits;
mod util;

#[cfg(test)]
mod testing;

pub use counts::*;

/// Type alias for a [`Result`] where the
/// error type is [`PopgenError`]
pub type PopgenResult<T> = Result<T, PopgenError>;

#[non_exhaustive]
#[derive(Debug)]
/// Error type
pub enum PopgenError {
    /// Wraps [`std::io::Error`].
    Io(std::io::Error),
    /// Returned when [`Count`] values are invalid.
    NegativeCount(Count),
    /// When the sum of allele counts at a site is >
    /// than the sample size (in genomes) at the same site.
    TotalAllelesDeficient,
    /// When non-empty site count data structure are required.
    EmptySiteCounts,
    /// When two objects selected for combination with [`TryReduce`](crate::traits::TryReduce) do not have matching lengths along some dimension.
    MismatchedLengths(usize, usize),
    /// Catch-all for numeric errors
    CalculationError,
    /// Invalid sample set identifier
    InvalidSampleSet,
    /// General library error containing a string-like message.
    LibraryError(String),
}

impl std::fmt::Display for PopgenError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            PopgenError::NegativeCount(c) => {
                write!(f, "inputted allele count may not be negative; got {}", c)
            }
            PopgenError::TotalAllelesDeficient => write!(
                f,
                "stated total alleles is less than sum of counts of present variants"
            ),
            PopgenError::Io(e) => write!(f, "io error: {}", e),
            PopgenError::MismatchedLengths(a, b) => {
                write!(
                    f,
                    "structures were expected to be of the same length; {a} != {b}"
                )
            }
            PopgenError::EmptySiteCounts => write!(f, "empty site count data"),
            PopgenError::CalculationError => write!(f, "calculation produced an invalid value"),
            PopgenError::InvalidSampleSet => write!(f, "invalid sample set label or index"),
            PopgenError::LibraryError(msg) => write!(f, "{msg}"),
        }
    }
}

impl std::error::Error for PopgenError {}

impl From<std::io::Error> for PopgenError {
    fn from(e: std::io::Error) -> Self {
        PopgenError::Io(e)
    }
}

/// The index/identifier of an allele at a given site.
/// This type is primarily used when reading input data,
/// such as allele records from a VCF file, etc..
/*
This could be u32, but it isn't for these reasons:
- The upside would be an increase in memory bandwidth and the ability to fit more counts data (double the amount on a 64-bit system) in a cache line.
- The downside is that ergonomics become rougher.
  Indexing an array by AlleleID would require .into() to extract the number, then casting to usize, which is fallible, because it may narrow to 16 bits.
- The upside is not useful until data is quite large, so we lean to use usize, for now.
*/
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct AlleleID(usize);

impl From<usize> for AlleleID {
    fn from(value: usize) -> Self {
        Self(value)
    }
}

impl From<AlleleID> for usize {
    fn from(value: AlleleID) -> Self {
        value.0
    }
}

impl AlleleID {
    /// Get the inner value of this ID.
    pub fn inner(&self) -> usize {
        self.0
    }
}

/// Get the crate version
pub fn version() -> &'static str {
    env!("CARGO_PKG_VERSION")
}
