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

fn pairwise_diffs(alleles: &[usize]) -> (i64, i64) {
    let mut num_differences = 0_i64;
    let mut num_comparisons = 0_i64;
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
    (num_differences, num_comparisons)
}

fn pi_site(genotypes: &mut dyn Iterator<Item = GenotypeData>) -> f64 {
    let alleles = flatten_to_alleles(genotypes);
    let (num_differences, num_comparisons) = pairwise_diffs(&alleles);
    num_differences as f64 / num_comparisons as f64
}

// O(N^2) implementation of the Nei/Tajima diversity measure.
pub fn pi<'s>(sites: &'s mut dyn Iterator<Item = &'s mut Site>) -> f64 {
    sites.map(|s| pi_site(&mut s.iter().cloned())).sum::<f64>()
}

fn watterson_theta_denominator(n: usize) -> f64 {
    let mut rv = 0.;
    for i in 1..n {
        rv += 1. / (i as f64)
    }
    rv
}

pub fn watterson_theta<'s>(sites: &'s mut dyn Iterator<Item = &'s mut Site>) -> f64 {
    let mut rv = 0.0;
    for s in sites {
        let alleles = flatten_to_alleles(&mut s.iter().cloned());
        let an = watterson_theta_denominator(alleles.len());
        let num_unique_alleles = alleles
            .into_iter()
            .collect::<std::collections::HashSet<_>>()
            .len();
        if num_unique_alleles > 1 {
            rv += (num_unique_alleles as f64 - 1.) / an;
        }
    }
    rv
}
