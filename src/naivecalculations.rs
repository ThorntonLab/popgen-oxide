use crate::testdata::GenotypeData;
use crate::testdata::Site;
use crate::MultiSiteCounts;

pub fn pi_site(genotypes: &mut dyn Iterator<Item = GenotypeData>) -> f64 {
    todo!()
}

pub fn pi<'s>(sites: &'s mut dyn Iterator<Item = &'s mut Site>) -> f64 {
    sites.map(|s| pi_site(&mut s.iter().cloned())).sum::<f64>()
}
