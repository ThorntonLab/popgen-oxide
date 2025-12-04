use std::fmt::Debug;
use std::str::FromStr;

pub mod adapter;
mod counts;
#[cfg(feature = "tskit")]
mod from_tree_sequence;
pub mod iter;
pub mod stats;
mod util;

#[cfg(test)]
mod naivecalculations;
#[cfg(test)]
mod test;
#[cfg(test)]
mod testdata;

pub use counts::*;

pub type PopgenResult<T> = Result<T, PopgenError>;

#[cfg(feature = "tskit")]
pub use tskit;

#[cfg(feature = "tskit")]
pub use from_tree_sequence::FromTreeSequenceOptions;

#[non_exhaustive]
#[derive(Debug)]
pub enum PopgenError {
    #[cfg(feature = "noodles")]
    NoodlesVCF(std::io::Error),
    #[cfg(feature = "tskit")]
    Tskit(tskit::TskitError),
    Io(std::io::Error),
    NegativeCount(Count),
    TotalAllelesDeficient,
    MismatchedSliceLength,
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
            #[cfg(feature = "tskit")]
            PopgenError::Tskit(e) => write!(f, "tskit error: {}", e),
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
