use std::fmt::Debug;
use std::str::FromStr;

pub mod adapter;
pub mod counts;
#[cfg(feature = "tskit")]
mod from_tree_sequence;
pub mod iter;
pub mod stats;
mod test;
pub(crate) mod util;

pub type PopgenResult<T> = Result<T, PopgenError>;

#[cfg(feature = "tskit")]
pub use tskit;

#[cfg(feature = "tskit")]
pub use from_tree_sequence::FromTreeSequenceOptions;

#[non_exhaustive]
#[derive(Debug, thiserror::Error)]
pub enum PopgenError {
    #[cfg(feature = "noodles")]
    #[error("couldn't handle VCF: {0}")]
    NoodlesVCF(#[from] std::io::Error),
    #[cfg(feature = "tskit")]
    #[error("tskit error: {0}")]
    Tskit(#[from] tskit::TskitError),
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

type Count = i64;
