use crate::iter::SiteCounts;
use crate::Count;
use std::iter::Sum;

/// A statistic calculable from and applicable to one site/locus.
pub trait SiteStatistic {
    fn from_site(site: SiteCounts) -> Self;
}

/// A statistic calculable from and applicable to a collection of sites/loci.
pub trait GlobalStatistic {
    fn from_iter_sites<'counts, I>(iter: I) -> Self
    where
        I: Iterator<Item = SiteCounts<'counts>>;
}

/// The chance pi that two uniformly randomly chosen genotypes at a site are different.
/// Expect a value in the 1e-3 range.
/// The complement of this probability is the site homozygosity.
///
/// This is known roughly as site heterozygosity or site diversity.
#[derive(Debug)]
pub struct Pi(pub f64);

impl Sum for Pi {
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        Self(iter.map(|x| x.0).sum())
    }
}

impl SiteStatistic for Pi {
    fn from_site(site: SiteCounts) -> Self {
        // technically should divide both by two here and below but it cancels out
        let num_pairs = {
            let count: i64 = site.counts.iter().sum();
            count * (count - 1)
        };

        // the number of pairs where the two samples are homozygous, summed over every genotype
        let num_homozygous_pairs: Count = site.counts.iter().map(|count| count * (count - 1)).sum();

        Self(1f64 - (num_homozygous_pairs as f64 / num_pairs as f64))
    }
}

/// The expected number of differences between two samples over all sites, the "expected pairwise diversity".
///
/// This is the sum of [`Pi`] over all sites.
#[derive(Debug)]
pub struct GlobalPi(pub f64);

impl GlobalStatistic for GlobalPi {
    fn from_iter_sites<'counts, I>(iter: I) -> Self
    where
        I: Iterator<Item = SiteCounts<'counts>>,
    {
        Self(iter.map(Pi::from_site).sum::<Pi>().0)
    }
}