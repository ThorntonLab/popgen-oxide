use crate::iter::SiteCounts;
use crate::Count;
use std::iter::Sum;

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
        self.0 += Pi::from_site(site).0;
    }

    fn as_raw(&self) -> f64 {
        self.0
    }
}