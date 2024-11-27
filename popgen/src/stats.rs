use std::cmp::max;
use crate::iter::SiteCounts;
use crate::Count;

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
        while let Some(w) = iter.next() {
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
