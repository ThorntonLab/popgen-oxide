use crate::util::StrictlyLowerTriangular;
use crate::{
    AlleleCounts, Count, MultiSampleAlleleCounts, PopgenError, PopgenResult, SampleAlleleCounts,
};
use std::cmp::max;

/// A statistic calculable from and applicable to one site/locus.
pub trait SiteStatistic {
    fn from_site(site: AlleleCounts) -> Self;
    fn as_raw(&self) -> f64;
}

/// A statistic calculable from and applicable to a collection of sites/loci.
pub trait GlobalStatistic {
    /// Instantiate a `Self` from an iterator over [`AlleleCounts`].
    ///
    /// The implementation depends on [`GlobalStatistic::try_add_site`].
    ///
    /// # Errors
    ///
    /// * If the iterator is empty, implementors may return [`PopgenError::EmptySiteCounts`],
    ///   because many statistics are not meaningfully defined over 0 sites.
    /// * Any error resulting from adding a single site ([`Self::try_add_site`]) may also be raised.
    fn try_from_iter_sites<'counts, I>(iter: I) -> Result<Self, PopgenError>
    where
        I: Iterator<Item = AlleleCounts<'counts>>,
        Self: Default,
    {
        let mut p = iter.peekable();
        if p.peek().is_some() {
            let mut ret = Self::default();
            for site in p {
                ret.try_add_site(site)?
            }
            Ok(ret)
        } else {
            Err(PopgenError::EmptySiteCounts)
        }
    }

    /// Update the value of a statistic from a [`AlleleCounts`]
    ///
    /// # Errors
    ///
    /// * Implementations should return [`PopgenError::CalculationError`]
    ///   if updating results in an invalid value.
    ///
    /// # Notes
    ///
    /// In general, one cannot assume that an update can be rolled
    /// back in the event of an error.
    ///
    /// [`SampleAlleleCounts`] and [`crate::MultiSampleAlleleCounts`] are designed
    /// such that `site` cannot contain empty data. However, it is valuable
    /// for implementations of this function to at least do the following:
    /// ```no_compile
    /// debug_assert!(!site.counts().is_empty());
    /// ```
    fn try_add_site(&mut self, site: AlleleCounts) -> Result<(), PopgenError>;
    fn as_raw(&self) -> f64;
}

/// A [`GlobalStatistic`], with the following additional guarantees:
/// - The statistic can be updated with at least one site when given a reference to another `Self` to which site(s) have already been added.
/// - This composition remains correct under reordering (commutation and association).
///
/// It is also possible to implement this trait if a sufficient *portion* of computation can be done
/// under the above conditions.
/// Such types can use [`GlobalStatistic::as_raw`] to perform inexpensive finalizing computations
/// if needed.
/// For an example of this use case, see [`TajimasD`].
///
/// # Errors
///
/// - The statistic is not required to be in a valid state after [`try_combine`](Self::try_combine) fails.
///
/// # Example
/// ```
/// # use popgen::{PopgenResult, AlleleCounts};
/// # use popgen::stats::{GlobalStatistic, SiteComposable};
/// #
/// // a very simple statistic
/// #[derive(Default)]  // to define the result over 0 sites
/// struct NumSites(u64);
///
/// impl SiteComposable for NumSites {
///     fn try_combine(&mut self, other: &Self) -> PopgenResult<()> {
///         self.0 += other.0;
///         Ok(())
///     }
/// }
///
/// impl GlobalStatistic for NumSites {
///     fn try_add_site(&mut self, site: AlleleCounts) -> PopgenResult<()> {
///         self.0 += 1;
///         Ok(())
///     }
///
///     fn as_raw(&self) -> f64 {
///         self.0 as f64
///     }
/// }
/// ```
pub trait SiteComposable: GlobalStatistic {
    fn try_combine(&mut self, other: &Self) -> crate::PopgenResult<()>;
}

pub fn windowed<Stat, GetWindow, E>(
    window_size: i64,
    stride: i64,
    start_pos: i64,
    end_pos: i64,
    mut get_window: GetWindow,
) -> PopgenResult<Result<Vec<Stat>, E>>
where
    Stat: SiteComposable + Default + Clone,
    GetWindow: FnMut(i64, i64) -> Result<SampleAlleleCounts, E>,
{
    let use_add_remove = window_size > stride;
    let mut intersection = None::<Stat>;

    let mut ret = Vec::with_capacity(((end_pos - start_pos + 1) / stride + 1) as usize);

    let mut window_start = start_pos;
    // TODO: parallelize this under some condition

    while window_start < end_pos {
        let mut accum = Stat::default();
        if !use_add_remove {
            let counts = match get_window(window_start, (window_start + window_size).min(end_pos)) {
                Ok(window) => window,
                Err(e) => return Ok(Err(e)),
            };
            for site in counts.iter() {
                accum.try_add_site(site)?;
            }
        } else {
            let new_intersection = match get_window(
                window_start + window_size - stride,
                window_start + window_size,
            ) {
                Ok(window) => window,
                Err(e) => return Ok(Err(e)),
            }
            .iter()
            .try_fold(Stat::default(), |mut stat, site| {
                stat.try_add_site(site)?;
                PopgenResult::Ok(stat)
            })?;

            if let Some(ref intersection) = intersection {
                accum.try_combine(intersection)?;

                let new_part = match get_window(window_start + stride, window_start + window_size) {
                    Ok(window) => window,
                    Err(e) => return Ok(Err(e)),
                };
                new_part
                    .iter()
                    .try_fold(Stat::default(), |mut stat, site| {
                        stat.try_add_site(site)?;
                        PopgenResult::Ok(stat)
                    })?;
            } else {
                let first_part = match get_window(window_start, window_start + window_size - stride)
                {
                    Ok(window) => window,
                    Err(e) => return Ok(Err(e)),
                };
                first_part
                    .iter()
                    .try_fold(Stat::default(), |mut stat, site| {
                        stat.try_add_site(site)?;
                        PopgenResult::Ok(stat)
                    })?;

                accum.try_combine(&new_intersection)?;
                intersection = Some(new_intersection);
            }

            ret.push(accum);
            window_start += stride;
        }
    }

    Ok(Ok(ret))
}

/// The expected number of differences between two samples over all sites, the "expected pairwise diversity".
///
/// Note that this statistic is **not defined** over an empty dataset; the denominator is the number of valid pairwise comparisons, which is 0 in this case.
/// Users should only use the [`Default`] implementation if they plan to do updates after construction, or to deliberately replace [`PopgenError::EmptySiteCounts`].
#[derive(Debug, Copy, Clone, Default)]
#[repr(transparent)]
pub struct Diversity(f64);

impl GlobalStatistic for Diversity {
    fn try_add_site(&mut self, site: AlleleCounts) -> Result<(), PopgenError> {
        debug_assert!(!site.counts().is_empty());
        // technically should divide both by two here and below but it cancels out
        let num_pairs = {
            let count: i64 = site.counts.iter().sum();
            count * (count - 1)
        };

        // the number of pairs where the two samples are homozygous, summed over every genotype
        let num_homozygous_pairs: Count = site.counts.iter().map(|count| count * (count - 1)).sum();

        self.0 += 1f64 - (num_homozygous_pairs as f64 / num_pairs as f64);
        if self.0.is_nan() {
            Err(PopgenError::CalculationError)
        } else {
            Ok(())
        }
    }

    fn as_raw(&self) -> f64 {
        self.0
    }
}

impl SiteComposable for Diversity {
    fn try_combine(&mut self, other: &Self) -> PopgenResult<()> {
        if self.0.is_nan() || other.0.is_nan() {
            return Err(PopgenError::CalculationError);
        }
        self.0 += other.as_raw();
        Ok(())
    }
}

/// Watterson's theta: see [Watterson's article](https://doi.org/10.1016%2F0040-5809%2875%2990020-9) and [Wikipedia](https://en.wikipedia.org/wiki/Watterson_estimator)
///
/// This statistic will never return [`PopgenError::EmptySiteCounts`], because it has a well-defined meaning under 0 sites.
/// The value returned by [`Default`] carries this meaning.
#[derive(Debug, Copy, Clone, Default)]
#[repr(transparent)]
pub struct WattersonTheta(f64);

impl GlobalStatistic for WattersonTheta {
    fn try_add_site(&mut self, site: AlleleCounts) -> Result<(), PopgenError> {
        debug_assert!(!site.counts().is_empty());
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

        if num_variants != 1 {
            let harmonic = (1..total_samples).map(|i| 1f64 / i as f64).sum::<f64>();
            self.0 += (num_variants - 1) as f64 / harmonic;
        }

        Ok(())
    }

    fn as_raw(&self) -> f64 {
        self.0
    }
}

impl SiteComposable for WattersonTheta {
    fn try_combine(&mut self, other: &Self) -> PopgenResult<()> {
        self.0 += other.0;
        Ok(())
    }
}

/// Tajima's D, as proposed in [Tajima 1989](https://academic.oup.com/genetics/article/123/3/585/5998755?login=false).
/// See also [Wikipedia](https://en.wikipedia.org/wiki/Tajima%27s_D#Mathematical_details) for the equations restated.
///
/// Note that this statistic is **not defined** over an empty dataset, because it depends on [`Diversity`] (see that documentation).
/// Users should only use the [`Default`] implementation if they plan to do updates after construction, or to deliberately replace [`PopgenError::EmptySiteCounts`].
#[derive(Debug, Copy, Clone, Default)]
pub struct TajimasD {
    k_hat: Diversity,
    theta: WattersonTheta,
    num_samples: usize,
    num_sites: usize,
}

impl GlobalStatistic for TajimasD {
    fn try_add_site(&mut self, site: AlleleCounts) -> Result<(), PopgenError> {
        debug_assert!(!site.counts().is_empty());
        self.k_hat.try_add_site(site.clone())?;
        self.theta.try_add_site(site.clone())?;

        self.num_sites += 1;
        // this is not perfect but that's fine
        self.num_samples = max(self.num_samples, site.total_alleles as usize);
        Ok(())
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

impl SiteComposable for TajimasD
where
    Diversity: SiteComposable,
    WattersonTheta: SiteComposable,
{
    fn try_combine(&mut self, other: &Self) -> PopgenResult<()> {
        self.k_hat.try_combine(&other.k_hat)?;
        self.theta.try_combine(&other.theta)?;
        self.num_samples += other.num_samples;
        self.num_sites += other.num_sites;
        Ok(())
    }
}

/// Fixation statistics as in [Charlesworth (1998)](https://doi.org/10.1093/oxfordjournals.molbev.a025953) and [Peter, 2016](https://pubmed.ncbi.nlm.nih.gov/26857625/).
///
/// Construction of this type from an arbitrary collection of [`SampleAlleleCounts`] is not sound,
/// because the invariant of [`crate::counts::MultiSampleAlleleCounts`] is required.
#[allow(non_camel_case_types)]
#[allow(non_snake_case)]
#[derive(Clone, Debug)]
pub struct FStatistics {
    /// (population number, weight) pairs
    populations: Vec<(usize, f64)>,
    // total pi_t derivable from other terms, no need to store anything new
    /// for pi_s
    diversity_within: Vec<f64>,
    // keep numerator and denominator apart for incremental update
    pi_s: (f64, f64),
    /// for pi_b
    divergence_between: StrictlyLowerTriangular<f64>,
    pi_b: (f64, f64),
}

impl FStatistics {
    fn new() -> Self {
        Self {
            populations: vec![],
            diversity_within: vec![],
            pi_s: (0.0, 0.0),
            divergence_between: StrictlyLowerTriangular::new(),
            pi_b: (0.0, 0.0),
        }
    }

    /// Add a population and its weight for this statistic.
    /// It is assumed that the inputted weight(s) sum to 1.
    ///
    /// # Errors
    /// See [`crate::stats::GlobalPi`].
    fn try_add_population(
        &mut self,
        populations: &MultiSampleAlleleCounts,
        population_num: usize,
        weight: f64,
    ) -> Result<(), PopgenError> {
        let diversity_new_site =
            Diversity::try_from_iter_sites(populations.iter_sites_in(population_num))?.as_raw();
        self.diversity_within.push(diversity_new_site);

        self.pi_s.0 += weight * weight * diversity_new_site;
        self.pi_s.1 += weight * weight;

        // there are more possible pairs of populations now
        self.divergence_between
            .try_extend(
                self.populations
                    .iter()
                    .map(|(existing_pop, existing_pop_weight)| {
                        let divergence_ij = populations
                            .iter_sites_in(*existing_pop)
                            .zip(populations.iter_sites_in(population_num))
                            .map(|(s1, s2)| {
                                if s1.total_alleles == 0 || s2.total_alleles == 0 {
                                    return Err(PopgenError::EmptySiteCounts);
                                }

                                // do complement of diversity, i.e. expected homozygosity

                                let total_comparisons = (s1.counts().iter().sum::<Count>()
                                    * s2.counts().iter().sum::<Count>())
                                    as i32;
                                if total_comparisons == 0 {
                                    return Err(PopgenError::EmptySiteCounts);
                                }

                                let num_homozygous = (0..max(s1.counts.len(), s2.counts.len()))
                                    .map(|variant_num| {
                                        // how many homozygous pairs?
                                        s1.counts.get(variant_num).unwrap_or(&0)
                                            * s2.counts.get(variant_num).unwrap_or(&0)
                                    })
                                    .sum::<i64>();

                                Ok(1. - num_homozygous as f64 / (total_comparisons as f64))
                            })
                            .sum::<Result<f64, PopgenError>>()?;

                        self.pi_b.0 += weight * existing_pop_weight * divergence_ij;
                        self.pi_b.1 += weight * existing_pop_weight;

                        PopgenResult::Ok(divergence_ij)
                    }),
            )?;

        self.populations.push((population_num, weight));
        Ok(())
    }

    /// Stream selected populations of a [`MultiSampleAlleleCounts`] into a computation of [`FStatistics`].
    ///
    /// Populations are both selected for inclusion/exclusion and assigned a weight using the input `pred`, which is called with the index of a population.
    /// The newly created struct immutably borrows from `self`.
    /// # Errors
    /// - If no populations are selected for inclusion.
    /// - If any population selected for inclusion has no sites or if any site on that population has zero present or total alleles.
    pub fn try_from_populations(
        populations: &MultiSampleAlleleCounts,
        mut pred: impl FnMut(usize) -> Option<f64>,
    ) -> Result<Self, PopgenError> {
        let mut ret = Self::new();

        let mut any = false;
        for pop_i in 0..populations.num_populations() {
            if let Some(weight) = pred(pop_i) {
                ret.try_add_population(populations, pop_i, weight)?;
                any = true;
            }
        }

        if !any {
            return Err(PopgenError::CalculationError);
        }

        Ok(ret)
    }

    ///.The two halves of the fraction [`Self::pi_s`].
    pub fn pi_s_parts(&self) -> (f64, f64) {
        self.pi_s
    }

    ///.The two halves of the fraction [`Self::pi_b`].
    pub fn pi_b_parts(&self) -> (f64, f64) {
        self.pi_b
    }

    /// The total diversity of these populations as defined by Charlesworth (1998) equations 1a and 2.
    pub fn pi_t(&self) -> f64 {
        let (pi_s_unweighted, _) = self.pi_s_parts();
        let (pi_b_unweighted, _) = self.pi_b_parts();
        pi_s_unweighted + 2. * pi_b_unweighted
    }

    /// The diversity of each population against itself as defined by Charlesworth (1998) equation 1b.
    /// [`None`] if no populations have been added so this fraction is undefined.
    pub fn pi_s(&self) -> Option<f64> {
        let pi_s = self.pi_s_parts();
        match pi_s.1 {
            0. => None,
            denom => Some(pi_s.0 / denom),
        }
    }

    /// The diversity between distinct populations as defined by Charlesworth (1998) equation 1c.
    /// [`None`] if no populations have been added so this fraction is undefined.
    pub fn pi_b(&self) -> Option<f64> {
        let pi_b = self.pi_b_parts();
        match pi_b.1 {
            0. => None,
            denom => Some(pi_b.0 / denom),
        }
    }

    /// As defined before equation 2.
    /// Calculated as pi_B - pi_S, from Charlesworth's pi_S + pi_D = pi_B.
    /// [`None`] if any of the required terms is undefined.
    pub fn pi_d(&self) -> Option<f64> {
        self.pi_b().zip(self.pi_s()).map(|(b, s)| b - s)
    }

    // pi_(T-S) is done fastest as pi_T - pi_S instead of with pi_D as in eqn 2b

    /// F_ST as defined by [Weir and Cockerham (1984)](https://doi.org/10.1111/j.1558-5646.1984.tb05657.x).
    /// [`None`] if any of the required terms is undefined.
    pub fn weir_cockerham(&self) -> Option<f64> {
        // eqn 3a
        self.pi_d().zip(self.pi_s()).map(|(d, s)| d / (s + d))
    }

    /// F_ST as defined by [Slatkin (1993)](https://doi.org/10.1111/j.1558-5646.1993.tb01215.x).
    /// [`None`] if any of the required terms is undefined.
    pub fn slatkin(&self) -> Option<f64> {
        // eqn 3b
        self.pi_d().zip(self.pi_s()).map(|(d, s)| d / (2. * s + d))
    }

    /// F_ST as defined by [Hudson, Boos, and Kaplan (1992)](https://doi.org/10.1093/oxfordjournals.molbev.a040703).
    /// [`None`] if any of the required terms is undefined.
    pub fn hudson_boos_kaplan(&self) -> Option<f64> {
        Some(self.pi_t()).zip(self.pi_s()).map(|(t, s)| (t - s) / t)
    }

    fn internal_index_for(&self, deme: usize) -> Option<usize> {
        match self.populations.len() {
            0..100 => self.populations.iter().position(|(p, _w)| p == &deme),
            _more => self
                .populations
                .binary_search_by_key(&deme, |(p, _w)| *p)
                .ok(),
        }
    }

    /// Get the diversity of this deme among its samples.
    /// The deme number must follow the indexes of demes used to create this type.
    ///
    /// # Errors
    ///
    //// * If `deme` is out of range, return [`PopgenError::InvalidDeme`]
    pub fn pi_within(&self, deme: usize) -> PopgenResult<f64> {
        Ok(self.diversity_within[self
            .internal_index_for(deme)
            .ok_or(PopgenError::InvalidDeme)?])
    }

    /// Get the [`Diversity`] of these two demes, comparing a sample from one against a sample from the other.
    /// The deme numbers must follow the indexes of demes used to create this type.
    ///
    /// # Errors
    ///
    /// * If `deme1` or `deme2` is out of range, return [`PopgenError::InvalidDeme`]
    pub fn pi_between(&self, deme1: usize, deme2: usize) -> PopgenResult<f64> {
        Ok(*self.divergence_between.get(
            self.internal_index_for(deme1)
                .ok_or(PopgenError::InvalidDeme)?,
            self.internal_index_for(deme2)
                .ok_or(PopgenError::InvalidDeme)?,
        ))
    }

    /// Calculate F2(deme1, deme2).
    /// The deme numbers must follow the indexes of demes used to create this type.
    ///
    /// We follow Equation 17 from
    /// [Peter, 2016](https://pubmed.ncbi.nlm.nih.gov/26857625/).
    ///
    /// # Errors
    ///
    /// * If `deme1` or `deme2` is out of range, return [`PopgenError::InvalidDeme`]
    pub fn f2(&self, deme1: usize, deme2: usize) -> Result<f64, PopgenError> {
        let deme1_internal = self
            .internal_index_for(deme1)
            .ok_or(PopgenError::InvalidDeme)?;
        let deme2_internal = self
            .internal_index_for(deme2)
            .ok_or(PopgenError::InvalidDeme)?;

        let divergence_12 = self.divergence_between.get(deme1_internal, deme2_internal);
        let diversity_11 = self
            .diversity_within
            .get(deme1_internal)
            .ok_or(PopgenError::InvalidDeme)?;
        let diversity_22 = self
            .diversity_within
            .get(deme2_internal)
            .ok_or(PopgenError::InvalidDeme)?;
        Ok(divergence_12 - (diversity_11 + diversity_22) / 2.)
    }

    /// Calculate F3(deme1; deme2, deme3).
    /// The deme numbers must follow the indexes of demes used to create this type.
    ///
    /// We follow Equation 20b from
    /// [Peter, 2016](https://pubmed.ncbi.nlm.nih.gov/26857625/),
    /// with some change of notation.
    /// He writes F3(deme X; deme 1, deme2).
    ///
    /// Originally due to Reich 2009 (as cited in Peter).
    ///
    /// # Errors
    ///
    /// * If any deme index is out of range, return [`PopgenError::InvalidDeme`]
    pub fn f3(&self, deme1: usize, deme2: usize, deme3: usize) -> Result<f64, PopgenError> {
        let a = self.f2(deme1, deme2)?;
        let b = self.f2(deme1, deme3)?;
        let c = self.f2(deme2, deme3)?;
        Ok((a + b - c) / 2.)
    }

    /// Calculate F4(deme1, deme2; deme3, deme4).
    /// The deme numbers must follow the indexes of demes used to create this type.
    ///
    /// We follow Equation 24b from
    /// [Peter, 2016](https://pubmed.ncbi.nlm.nih.gov/26857625/).
    ///
    /// Originally due to Reich 2009 (as cited in Peter).
    ///
    /// # Errors
    ///
    /// * If any deme index is out of range, return [`PopgenError::InvalidDeme`]
    pub fn f4(
        &self,
        deme1: usize,
        deme2: usize,
        deme3: usize,
        deme4: usize,
    ) -> Result<f64, PopgenError> {
        let a = self.f2(deme1, deme4)?;
        let b = self.f2(deme2, deme3)?;
        let c = self.f2(deme1, deme3)?;
        let d = self.f2(deme2, deme4)?;
        Ok((a + b - c - d) / 2.)
    }
}
