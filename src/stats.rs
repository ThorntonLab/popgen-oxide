use crate::iter::SiteCounts;
use crate::util::UnorderedPair;
use crate::{Count, MultiSiteCounts};
use itertools::Itertools;
use std::cmp::max;
use std::collections::{HashMap, HashSet};

/// A statistic calculable from and applicable to one site/locus.
pub trait SiteStatistic {
    fn from_site(site: SiteCounts) -> Self;
    fn as_raw(&self) -> f64;
}

/// A statistic calculable from and applicable to a collection of sites/loci.
pub trait GlobalStatistic {
    fn from_iter_sites<'counts, I>(iter: I) -> Self
    where
        I: Iterator<Item = SiteCounts<'counts>>,
        Self: Default,
    {
        let mut ret = Self::default();
        for site in iter {
            ret.add_site(site)
        }
        ret
    }

    fn add_site(&mut self, site: SiteCounts);
    fn as_raw(&self) -> f64;
}

/// The expected number of differences between two samples over all sites, the "expected pairwise diversity".
///
/// This is the sum of [`Pi`] over all sites.
#[derive(Debug, Copy, Clone, Default)]
#[repr(transparent)]
pub struct GlobalPi(f64);

impl GlobalStatistic for GlobalPi {
    fn add_site(&mut self, site: SiteCounts) {
        // technically should divide both by two here and below but it cancels out
        let num_pairs = {
            let count: i64 = site.counts.iter().sum();
            count * (count - 1)
        };

        // the number of pairs where the two samples are homozygous, summed over every genotype
        let num_homozygous_pairs: Count = site.counts.iter().map(|count| count * (count - 1)).sum();

        self.0 += 1f64 - (num_homozygous_pairs as f64 / num_pairs as f64)
    }

    fn as_raw(&self) -> f64 {
        self.0
    }
}

impl From<&MultiSiteCounts> for GlobalPi {
    fn from(value: &MultiSiteCounts) -> Self {
        Self::from_iter_sites(value.iter())
    }
}

/// Watterson's theta: see [Watterson's article](https://doi.org/10.1016%2F0040-5809%2875%2990020-9) and [Wikipedia](https://en.wikipedia.org/wiki/Watterson_estimator)
#[derive(Debug, Copy, Clone, Default)]
#[repr(transparent)]
pub struct WattersonTheta(f64);

impl GlobalStatistic for WattersonTheta {
    fn add_site(&mut self, site: SiteCounts) {
        // trying our very hardest to encourage optimization and SIMD here
        // also optimizing with the typical two-element slice in mind
        let mut iter = site.counts.chunks_exact(2);
        let mut num_variants = 0;
        let mut total_samples = 0;
        for w in iter.by_ref() {
            // big idea: with chunks_exact this cast is infallible and zero-cost
            // this cast also enables use of 128-bit and SIMD instructions
            let w: &[i64; 2] = w.try_into().expect("slice with incorrect length");
            w.iter().filter(|&&c| c > 0).for_each(|&c| {
                num_variants += 1;
                total_samples += c;
            })
        }
        iter.remainder().iter().filter(|&&c| c > 0).for_each(|&c| {
            num_variants += 1;
            total_samples += c;
        });

        if num_variants > 1 {
            let harmonic = (1..total_samples).map(|i| 1f64 / i as f64).sum::<f64>();
            self.0 += (num_variants - 1) as f64 / harmonic;
        } else {
            // then this site isn't actually polymorphic; meh
        }
    }

    fn as_raw(&self) -> f64 {
        self.0
    }
}

/// Tajima's D, as proposed in [Tajima 1989](https://academic.oup.com/genetics/article/123/3/585/5998755?login=false).
/// See also [Wikipedia](https://en.wikipedia.org/wiki/Tajima%27s_D#Mathematical_details) for the equations restated.
#[derive(Debug, Copy, Clone, Default)]
pub struct TajimaD {
    k_hat: GlobalPi,
    theta: WattersonTheta,
    num_samples: usize,
    num_sites: usize,
}

impl GlobalStatistic for TajimaD {
    fn add_site(&mut self, site: SiteCounts) {
        self.k_hat.add_site(site.clone());
        self.theta.add_site(site.clone());

        self.num_sites += 1;
        // this is not perfect but that's fine
        self.num_samples = max(self.num_samples, site.total_alleles as usize);
    }

    fn as_raw(&self) -> f64 {
        // we are going to stick as closely as feasible to the exact nomenclature of the paper

        let n = self.num_samples;

        // eqn 3
        let a_1 = (1..n).map(|i| 1f64 / (i as f64)).sum::<f64>();
        // eqn 4
        let a_2 = (1..n).map(|i| 1f64 / ((i * i) as f64)).sum::<f64>();

        let n = n as f64;

        // eqn 8
        let b_1 = (n + 1.) / (3. * (n - 1.));
        // eqn 9
        let b_2 = (2. * (n * n + n + 3.)) / (9. * n * (n - 1.));

        // eqn 28
        let d = self.k_hat.as_raw() - self.theta.as_raw();

        // eqn 31
        let c_1 = b_1 - 1. / a_1;
        // eqn 32
        let c_2 = b_2 - (n + 2.) / (a_1 * n) + (a_2 / (a_1 * a_1));

        // eqn 36
        let e_1 = c_1 / a_1;
        // eqn 37
        let e_2 = c_2 / (a_1 * a_1 + a_2);

        #[allow(non_snake_case)]
        let S = self.num_sites as f64;

        // eqn 38
        #[allow(non_snake_case)]
        let D = d / (e_1 * S + e_2 * S * (S - 1.)).sqrt();
        D
    }
}

/// Fixation statistics as in [Charlesworth (1998)](https://doi.org/10.1093/oxfordjournals.molbev.a025953).
#[allow(non_camel_case_types)]
#[allow(non_snake_case)]
#[derive(Clone, Default, Debug)]
pub struct F_ST {
    /// (population, weight) pairs
    populations: Vec<(MultiSiteCounts, f64)>,
    // total pi_T derivable from other terms, no need to store anything new
    /// for pi_S
    diversity_within: Vec<f64>,
    // keep numerator and denominator apart for incremental update
    pi_S: (f64, f64),
    /// for pi_B
    diversity_between: HashMap<UnorderedPair<usize>, f64>,
    pi_B: (f64, f64),
}

// TODO: tests for all of these
impl F_ST {
    /// Construct a new instance of this statistic, ready to accept populations via [`Self::add_population`].
    pub fn new() -> Self {
        Self::default()
    }

    /// Add a population and its weight for this statistic.
    /// It is assumed that the inputted weight(s) sum to 1.
    pub fn add_population(&mut self, population: MultiSiteCounts, weight: f64) {
        let pi_new_site = GlobalPi::from(&population).as_raw();
        self.diversity_within.push(pi_new_site);

        self.pi_S.0 += weight * weight * pi_new_site;
        self.pi_S.1 += weight * weight;

        // there are more possible pairs of populations now
        for (i, (existing_pop, existing_site_weight)) in self.populations.iter().enumerate() {
            let pi_ij = existing_pop
                .iter()
                .zip(population.iter())
                .map(|(s1, s2)| {
                    1. -
                        // do complement of diversity, i.e. expected homozygosity
                        // for each variant...
                        (0..max(s1.counts.len(), s2.counts.len()))
                            .map(|variant_num| {
                                // how many homozygous pairs?
                                s1.counts.get(variant_num).unwrap_or(&0)
                                    * s2.counts.get(variant_num).unwrap_or(&0)
                            })
                            .sum::<i64>() as f64
                            // how many possible pairs?
                            / ((s1.total_alleles * s2.total_alleles) as f64)
                })
                .sum();

            // TODO: just use a linear vector?
            self.diversity_between
                .insert(UnorderedPair::new(i, self.populations.len()), pi_ij);

            self.pi_B.0 += weight * existing_site_weight * pi_ij;
            self.pi_B.1 += weight * existing_site_weight;
        }
        self.populations.push((population, weight));
    }

    /// Construct an [`F_STView`] from this struct with the provided `selected_populations`.
    /// It is assumed these are valid indices into the populations already provided to `self`, panicking otherwise.
    /// The terms necessary to compute F-statistics are recomputed immediately.
    pub fn only_populations(&self, selected_populations: HashSet<usize>) -> F_STView {
        #[allow(non_snake_case)]
        let mut pi_S = (0., 0.);
        #[allow(non_snake_case)]
        let mut pi_B = (0., 0.);

        for index in selected_populations.iter() {
            let (_, weight) = self.populations.get(*index).unwrap();
            pi_S.0 += weight * weight * self.diversity_within[*index];
            pi_S.1 += weight * weight;
        }

        for pair in selected_populations
            .iter()
            .tuple_combinations()
            .map(|(&i, &j)| UnorderedPair::new(i, j))
            .unique()
        {
            let UnorderedPair(i, j) = pair;
            pi_B.0 += self.populations[i].1
                * self.populations[j].1
                * self.diversity_between.get(&pair).unwrap();
            pi_B.1 += self.populations[i].1 * self.populations[j].1;
        }

        F_STView {
            inner: self,
            selected_populations,
            pi_S,
            pi_B,
        }
    }
}

pub trait FStatisticParts {
    #[allow(non_snake_case)]
    fn pi_S_parts(&self) -> (f64, f64);

    #[allow(non_snake_case)]
    fn pi_B_parts(&self) -> (f64, f64);

    /// The total diversity of these populations as defined by equations 1a and 2.
    #[allow(non_snake_case)]
    fn pi_T(&self) -> f64 {
        let (pi_S_unweighted, _) = self.pi_S_parts();
        let (pi_B_unweighted, _) = self.pi_B_parts();
        pi_S_unweighted + 2. * pi_B_unweighted
    }

    /// The diversity of each population against itself as defined by equation 1b.
    /// [`None`] if no populations have been added so this fraction is undefined.
    #[allow(non_snake_case)]
    fn pi_S(&self) -> Option<f64> {
        let pi_S = self.pi_S_parts();
        match pi_S.1 {
            0. => None,
            denom => Some(pi_S.0 / denom),
        }
    }

    /// The diversity between distinct populations as defined by equation 1c.
    /// [`None`] if no populations have been added so this fraction is undefined.
    #[allow(non_snake_case)]
    fn pi_B(&self) -> Option<f64> {
        let pi_B = self.pi_B_parts();
        match pi_B.1 {
            0. => None,
            denom => Some(pi_B.0 / denom),
        }
    }

    /// As defined before equation 2.
    /// Calculated as pi_B - pi_S, from Charlesworth's pi_S + pi_D = pi_B.
    /// [`None`] if any of the required terms is undefined.
    #[allow(non_snake_case)]
    fn pi_D(&self) -> Option<f64> {
        self.pi_B().zip(self.pi_S()).map(|(b, s)| b - s)
    }

    // pi_(T-S) is done fastest as pi_T - pi_S instead of with pi_D as in eqn 2b
}

pub trait FStatistics: FStatisticParts {
    /// F_ST as defined by [Weir and Cockerham (1984)](https://doi.org/10.1111/j.1558-5646.1984.tb05657.x).
    /// [`None`] if any of the required terms is undefined.
    fn weir_cockerham(&self) -> Option<f64> {
        // eqn 3a
        self.pi_D().zip(self.pi_S()).map(|(d, s)| d / (s + d))
    }

    /// F_ST as defined by [Slatkin (1993)](https://doi.org/10.1111/j.1558-5646.1993.tb01215.x).
    /// [`None`] if any of the required terms is undefined.
    fn slatkin(&self) -> Option<f64> {
        // eqn 3b
        self.pi_D().zip(self.pi_S()).map(|(d, s)| d / (2. * s + d))
    }

    /// F_ST as defined by [Hudson, Boos, and Kaplan (1992)](https://doi.org/10.1093/oxfordjournals.molbev.a040703).
    /// [`None`] if any of the required terms is undefined.
    fn hudson_boos_kaplan(&self) -> Option<f64> {
        Some(self.pi_T()).zip(self.pi_S()).map(|(t, s)| (t - s) / t)
    }
}

impl FStatisticParts for F_ST {
    #[allow(non_snake_case)]
    fn pi_S_parts(&self) -> (f64, f64) {
        self.pi_S
    }

    #[allow(non_snake_case)]
    fn pi_B_parts(&self) -> (f64, f64) {
        self.pi_B
    }
}

/// A view into an underlying [`F_ST`] with a subset of populations from which to compute F-statistics.
/// Constructed via [`F_ST::only_populations`].
#[allow(non_camel_case_types, non_snake_case)]
#[derive(Clone, Debug)]
pub struct F_STView<'a> {
    inner: &'a F_ST,
    selected_populations: HashSet<usize>,
    pi_S: (f64, f64),
    pi_B: (f64, f64),
}

impl<'a> F_STView<'a> {
    pub fn inner(&self) -> &'a F_ST {
        self.inner
    }

    pub fn selected_populations(&self) -> &HashSet<usize> {
        &self.selected_populations
    }
}

impl FStatisticParts for F_STView<'_> {
    #[allow(non_snake_case)]
    fn pi_S_parts(&self) -> (f64, f64) {
        self.pi_S
    }

    #[allow(non_snake_case)]
    fn pi_B_parts(&self) -> (f64, f64) {
        self.pi_B
    }
}
