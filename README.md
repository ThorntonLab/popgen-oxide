# Efficient types and methods for population genetics

This crate provides data structures and interfaces for efficient calculation from genetic variation data.

## Example

```rust
// Interface for a statistic that aggregates over a function
// of allele counts across sites with no need to distinguish
// ancestral/derived states.
use popgen::stats::UnpolarisedSiteStat;

// Interface to obtain the representation of 
// a statistic as a low-level type, providing
// the `as_raw` function used below.
use popgen::stats::StatRepresentation;

fn main() {
    // Raw variation data at two sites.
    // There are two alleles at the first site -- 
    // one copy of allele 0 and one copy of allele 2.
    // At the second site we have one copy each of alleles
    // 0 and 1.
    // (A None value would represent missing data.
    //  There are NO MAGIC NUMBERS like -1 used!)
    let data: Vec<Vec<Option<popgen::AlleleID>>> = vec![
        vec![Some(0.into()), Some(2.into())],
        vec![Some(0.into()), Some(1.into())],
    ];
    // Build our data structure representing allele counts
    // over all the sites.
    let counts = popgen::SampleAlleleCounts::try_from_tabular(data).unwrap();
    // Calculate diversity, aka "pi", aka "mean number of pairwise differences"
    let diversity =
        popgen::stats::Diversity::try_from_iter_sites(counts.iter_sample_set(0).unwrap()).unwrap();
    assert_eq!(diversity.as_raw(), 2.0)
}
```
