use crate::testdata::GenotypeData;
use crate::testdata::Site;
use crate::MultiSiteCounts;

pub fn pi_site(genotypes: &dyn Iterator<Item = GenotypeData>) -> f64 {
    todo!()
}

pub fn pi<'s>(sites: &'s dyn Iterator<Item = &'s Site>) -> f64 {
    todo!()
}

pub fn single_pop_counts<'s>(sites: &'s dyn Iterator<Item = &'s Site>) -> MultiSiteCounts {
    todo!()
}
