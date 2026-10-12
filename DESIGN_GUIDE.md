# Design guide

## General considerations

### Avoid "syntax sugar"

We want to avoid public API elements that are simple wrappers around existing functionality.
For example, if we have a function whose return type is `impl Iterator<Item=usize>`, then we do **not** want a second function that returns the same output as a `Vec<usize>`.
The reason is that the second function would simply be a dispatch to `Iterator::collect`, and is thus a kind of "convenience" function.
Such convenience functions are not without cost!
They have to be compiled, etc., during CI.
They are continuously compiled by LSPs during development of the code base.
They require documentation and maintenance.
Etc..

## Summary statistics

### Single locus/site statistics

Statistics must:

* be newtypes or structs that implement the relevant API traits for calculation from data

Statistics should:

* implement Debug, Copy and/or Clone as appropriate, and Default if appropriate.
  The derive macro implementations should be used if acceptable.
* implement Display, Eq/Ord or PartialEq/PartialOrd as appropriate for their underlying representation.
* implement From<Self> for conversion into their underlying representation

Some statistics require composititions of multiple values to calculate a final statistic.
An example is Tajima's D, which relies of diversity and Watterson's estimator of theta.
Another example are F statistics like Fst, F2, etc..

For these statistics, we want:

* A "builder" type that can implement the API traits for calculation.
  Example: FStatistics
* The builder type should have an API to access the underlying components.
  Example: FStatistics gives access to the stored diversity and divergence values.
* The builder type should have an API to resolve the trait into its final values.
  Example: TajimasDBuilder to Tajima's D via From.


