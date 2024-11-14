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
pub struct WattersonsTheta(f64);

impl GlobalStatistic for WattersonsTheta {
    fn add_site(&mut self, site: SiteCounts) {
        // trying our very hardest to encourage optimization and SIMD here
        // also optimizing with the typical two-element slice in mind
        let mut iter = site.counts.chunks_exact(2);
        let mut num_sites = 0;
        let mut total_samples = 0;
        while let Some(w) = iter.next() {
            // big idea: with chunks_exact this cast is infallible and zero-cost
            // this cast also enables use of 128-bit and SIMD instructions
            let w: &[i64; 2] = w.try_into().expect("slice with incorrect length");
            w.iter().filter(|&&c| c > 0).for_each(|&c| {
                num_sites += 1;
                total_samples += c;
            })
        }
        iter.remainder().iter().filter(|&&c| c > 0).for_each(|&c| {
            num_sites += 1;
            total_samples += c;
        });

        if num_sites > 1 {
            let harmonic = (1..total_samples).map(|i| 1f64 / i as f64).sum::<f64>();
            self.0 += (num_sites - 1) as f64 / harmonic;
        } else {
            // then this site isn't actually polymorphic; meh
        }
    }

    fn as_raw(&self) -> f64 {
        self.0
    }
}
