# Allele count data and summary statistics

## Example

```rust
// Import the trait
use popgen::stats::UnpolarisedSiteStat;
// Create some allele count data:
// Allele 0 occurs 10 times.
// Allele 1 occurs twice.
let c = vec![10_i64, 2];
// Create an instance of AlleleCounts by borrowing
// the above data.
// The second argument is the total number of alleles in
// the experiment.
// Here 14 is larger than the sum of elements in "c", which
// means that there were two missing allele calls in the input data
// set.
let counts = popgen::AlleleCounts::new(&c, 14).unwrap();
// Calculate diversity ("pi", "sum of site heterozygosity").
// The function try_from_iter_sites is part of the UnpolarisedSiteStat trait API.
let d = popgen::stats::Diversity::try_from_iter_sites(std::iter::once(counts)).unwrap();
```
