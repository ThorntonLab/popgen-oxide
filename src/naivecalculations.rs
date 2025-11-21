use crate::testdata::GenotypeData;
use crate::testdata::Site;

fn flatten_to_alleles(genotypes: &mut dyn Iterator<Item = GenotypeData>) -> Vec<usize> {
    let mut alleles = vec![];
    for gi in genotypes {
        for allele in gi.iter().flatten() {
            alleles.push(allele);
        }
    }
    alleles
}

fn pi_site(genotypes: &mut dyn Iterator<Item = GenotypeData>) -> f64 {
    let mut num_differences = 0_i64;
    let mut num_comparisons = 0_i64;
    let alleles = flatten_to_alleles(genotypes);
    for (i, j) in alleles.iter().enumerate() {
        for k in alleles.iter().skip(i + 1) {
            if j != k {
                num_differences += 1;
            }
            num_comparisons += 1;
        }
    }
    if !alleles.is_empty() {
        assert_eq!(
            num_comparisons as usize,
            alleles.len() * (alleles.len() - 1) / 2
        );
    }
    num_differences as f64 / num_comparisons as f64
}

// O(N^2) implementation of the Nei/Tajima diversity measure.
pub fn pi<'s>(sites: &'s mut dyn Iterator<Item = &'s mut Site>) -> f64 {
    sites.map(|s| pi_site(&mut s.iter().cloned())).sum::<f64>()
}
