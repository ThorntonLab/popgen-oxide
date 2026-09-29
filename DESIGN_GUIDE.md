# Design guide

## Summary statistics

### Single locus/site statistics

Statistics must:

* be newtypes or structs that implement the relevant API traits for calculation from data

Statistics should:

* implement Debug, Copy and/or Clone as appropriate, and Default if appropriate.
  The derive macro implementations should be used if acceptable.
* implement Display, Eq/Ord or PartialEq/PartialOrd as appropriate for their underlying representation.
* implement From<Self> for their underlying representation

Some statistics require composititions of multiple values to calculate a final statistic.
An examaple is Tajima's D, which relies of diversity and Watterson's estimator of theta.
Another example are F statistics like Fst, F2, etc..

For these statistics, we want:

* A "builder" type that can implement the API traits for calculation.
  Example: FStatistics
* The builder type should have an API to access the underlying components.
  Example: FStatistics gives access to the stored diversity and divergence values.
* The builder type should have an API to resolve the trait into its final values.
  Example: TajimasDBuilder to Tajima's D via From.


