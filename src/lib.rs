#![deny(rustdoc::broken_intra_doc_links)]
#![cfg_attr(doc_cfg, feature(doc_cfg))]

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

pub type VaristatResult<T> = Result<T, VaristatError>;

#[cfg(feature = "tskit")]
pub mod from_tskit {
    pub use super::from_tree_sequence::FromTreeSequenceError;
    pub use super::from_tree_sequence::FromTreeSequenceOptions;
}

#[non_exhaustive]
#[derive(Debug)]
pub enum VaristatError {
    #[cfg(feature = "noodles")]
    NoodlesVCF(std::io::Error),
    #[cfg(feature = "tskit")]
    Tskit(crate::from_tskit::FromTreeSequenceError),
    Io(std::io::Error),
    NegativeCount(Count),
    TotalAllelesDeficient,
    MismatchedSliceLength,
    EmptySiteCounts,
    CalculationError,
    InvalidDeme,
    LibraryError(String),
}

impl std::fmt::Display for VaristatError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            VaristatError::NegativeCount(c) => {
                write!(f, "inputted allele count may not be negative; got {}", c)
            }
            VaristatError::TotalAllelesDeficient => write!(
                f,
                "stated total alleles is less than sum of counts of present variants"
            ),
            VaristatError::Io(e) => write!(f, "io error: {}", e),
            VaristatError::MismatchedSliceLength => {
                write!(f, "slices were expected to be of the same length")
            }
            VaristatError::EmptySiteCounts => write!(f, "empty site count data"),
            VaristatError::CalculationError => write!(f, "calculation produced an invalid value"),
            VaristatError::InvalidDeme => write!(f, "invalid deme label or index"),
            VaristatError::LibraryError(msg) => write!(f, "{msg}"),
            #[cfg(feature = "tskit")]
            VaristatError::Tskit(e) => write!(f, "tskit conversion error: {}", e),
            #[cfg(feature = "noodles")]
            VaristatError::NoodlesVCF(e) => write!(f, "couldn't handle VCF: {}", e),
        }
    }
}

impl std::error::Error for VaristatError {}

impl From<std::io::Error> for VaristatError {
    fn from(e: std::io::Error) -> Self {
        VaristatError::Io(e)
    }
}

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
