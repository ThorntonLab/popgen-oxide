#![deny(rustdoc::broken_intra_doc_links)]
#![warn(missing_docs)]
#![cfg_attr(doc_cfg, feature(doc_cfg))]

//! Efficient types and interfaces for population genetics

use std::fmt::Debug;
use std::str::FromStr;

pub mod adapter;
mod counts;
#[cfg(feature = "tskit")]
mod from_tree_sequence;
pub mod iter;
pub mod stats;
pub mod traits;
mod util;

#[cfg(test)]
mod testing;

pub use counts::*;

/// Type alias for a [`Result`] where the
/// error type is [`PopgenError`]
pub type PopgenResult<T> = Result<T, PopgenError>;

#[cfg(feature = "tskit")]
#[cfg_attr(doc_cfg, doc(cfg(feature = "tskit")))]
pub mod from_tskit {
    pub use super::from_tree_sequence::FromTreeSequenceError;
    pub use super::from_tree_sequence::FromTreeSequenceOptions;
}

#[non_exhaustive]
#[derive(Debug)]
/// Error type
pub enum PopgenError {
    #[cfg_attr(doc_cfg, doc(cfg(feature = "noodles")))]
    #[cfg(feature = "noodles")]
    NoodlesVCF(std::io::Error),
    #[cfg(feature = "tskit")]
    #[cfg_attr(doc_cfg, doc(cfg(feature = "tskit")))]
    Tskit(crate::from_tskit::FromTreeSequenceError),
    Io(std::io::Error),
    /// Returned when [`Count`] values are invalid.
    NegativeCount(Count),
    TotalAllelesDeficient,
    MismatchedSliceLength,
    EmptySiteCounts,
    MismatchedSampleSetCount(usize, usize),
    CalculationError,
    InvalidDeme,
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
            PopgenError::MismatchedSliceLength => {
                write!(f, "slices were expected to be of the same length")
            }
            PopgenError::EmptySiteCounts => write!(f, "empty site count data"),
            PopgenError::MismatchedSampleSetCount(l, r) => write!(
                f,
                "cannot combine two collections with different sample set counts; {l} != {r}"
            ),
            PopgenError::CalculationError => write!(f, "calculation produced an invalid value"),
            PopgenError::InvalidDeme => write!(f, "invalid deme label or index"),
            PopgenError::LibraryError(msg) => write!(f, "{msg}"),
            #[cfg_attr(doc_cfg, doc(cfg(feature = "tskit")))]
            #[cfg(feature = "tskit")]
            PopgenError::Tskit(e) => write!(f, "tskit conversion error: {}", e),
            #[cfg_attr(doc_cfg, doc(cfg(feature = "noodles")))]
            #[cfg(feature = "noodles")]
            PopgenError::NoodlesVCF(e) => write!(f, "couldn't handle VCF: {}", e),
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

/// Get the crate version
pub fn version() -> &'static str {
    env!("CARGO_PKG_VERSION")
}
