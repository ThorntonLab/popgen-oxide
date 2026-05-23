mod naivecalculations;
mod testdata;

// tests of TYPES go below
mod multisitecounts;
#[cfg(feature = "tskit")]
mod test_try_from_tree_sequence;

// tests of CALCULATIONS go below

mod fst;
mod pi;
mod tajimasd;
mod wattersons_theta;

// tests of OTHER THINGS, like
// input formats, follow

mod concurrent;
#[cfg(feature = "noodles")]
mod noodles_vcf;
