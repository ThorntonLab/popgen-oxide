use std::collections::HashMap;

use crate::testing::testdata::GenotypeData;
use crate::testing::testdata::Site;

fn flatten_to_alleles(genotypes: &mut dyn Iterator<Item = GenotypeData>) -> Vec<usize> {
    let mut alleles = vec![];
    for gi in genotypes {
        for allele in gi.iter().flatten() {
            alleles.push(allele);
        }
    }
    alleles
}

fn flatten_to_allele_counts(
    genotypes: &mut dyn Iterator<Item = GenotypeData>,
) -> HashMap<usize, i64> {
    let alleles = flatten_to_alleles(genotypes);
    let mut counts: HashMap<usize, i64> = HashMap::new();
    for a in alleles {
        if let Some(c) = counts.get_mut(&a) {
            *c += 1
        } else {
            counts.insert(a, 1);
        }
    }
    counts
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
pub fn pi<'s>(sites: impl Iterator<Item = &'s Site>) -> f64 {
    let sites = sites.cloned();
    let mut sites = sites.peekable();
    if sites.peek().is_some() {
        sites.map(|s| pi_site(&mut s.iter().cloned())).sum::<f64>()
    } else {
        f64::NAN
    }
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

pub fn pi_ij(pop1: &mut dyn Iterator<Item = Site>, pop2: &mut dyn Iterator<Item = Site>) -> f64 {
    let mut accum = 0f64;

    for (site1, site2) in pop1.zip(pop2) {
        let mut homozygous = 0u64;
        let mut num_comparisons = 0u64;
        for gt1 in site1.iter().flat_map(|gt| gt.iter()) {
            for gt2 in site2.iter().flat_map(|gt| gt.iter()) {
                if let Some((l, r)) = gt1.zip(gt2) {
                    if l == r {
                        homozygous += 1;
                    }
                    num_comparisons += 1;
                }
            }
        }
        accum += 1.0 - (homozygous as f64 / num_comparisons as f64);
    }

    accum
}

// return pi_T, pi_S, pi_B
pub fn f_st<Sites>(populations: &mut dyn Iterator<Item = (f64, Sites)>) -> (f64, f64, f64)
where
    Sites: IntoIterator<Item = Site>,
{
    let pops_and_weights: Vec<(f64, Vec<Site>)> = populations
        .map(|(w, pop)| (w, pop.into_iter().collect()))
        .collect();
    let weights = pops_and_weights
        .iter()
        .map(|(w, _pop)| w)
        .copied()
        .collect::<Vec<_>>();
    let mut populations = pops_and_weights
        .into_iter()
        .map(|(_w, pop)| pop)
        .collect::<Vec<_>>();

    let pi_ii = populations
        .iter_mut()
        .map(|pop| pi(pop.iter()))
        .collect::<Vec<_>>();

    let mut pi_ij_computed = vec![vec![None; populations.len()]; populations.len()];
    for i in 0..populations.len() {
        for j in 0..i {
            pi_ij_computed[i][j] = Some(pi_ij(
                &mut populations[i].iter().cloned(),
                &mut populations[j].iter().cloned(),
            ));
            pi_ij_computed[j][i] = pi_ij_computed[i][j];
        }
    }

    // equation 1a
    let pi_total = {
        let pi_ii_term = (0..populations.len())
            .map(|i| weights[i] * weights[i] * pi_ii[i])
            .sum::<f64>();
        let pi_ij_term = {
            let mut tot = 0f64;
            for i in 0..populations.len() {
                for j in 0..i {
                    tot += weights[i] * weights[j] * pi_ij_computed[i][j].unwrap();
                }
            }
            2.0 * tot
        };

        pi_ii_term + pi_ij_term
    };

    // equation 1b
    let pi_self = ((0..populations.len())
        .map(|i| weights[i] * weights[i] * pi_ii[i])
        .sum::<f64>())
        / (0..populations.len())
            .map(|i| weights[i] * weights[i])
            .sum::<f64>();

    // equation 1c
    let pi_between = {
        let mut num = 0f64;
        let mut denom = 0f64;

        for i in 0..populations.len() {
            for j in 0..i {
                num += weights[i] * weights[j] * pi_ij_computed[i][j].unwrap();
                denom += weights[i] * weights[j];
            }
        }

        num / denom
    };

    (pi_total, pi_self, pi_between)
}

pub fn f2<Sites>(deme1: usize, deme2: usize, populations: &mut dyn Iterator<Item = Sites>) -> f64
where
    Sites: IntoIterator<Item = Site>,
{
    let mut naive = 0.0;
    let mut n = 0_usize;
    let freq_data: Vec<Vec<Site>> = populations.map(|pop| pop.into_iter().collect()).collect();
    for (d1, d2) in freq_data[deme1].iter().zip(freq_data[deme2].iter()) {
        let allele_counts1 = flatten_to_allele_counts(&mut d1.iter().cloned());
        let allele_counts2 = flatten_to_allele_counts(&mut d2.iter().cloned());
        let unique_allele_labels = {
            let mut temp = allele_counts1.keys().cloned().collect::<Vec<_>>();
            temp.extend(allele_counts2.keys().cloned());
            temp.sort_unstable();
            temp.dedup();
            temp
        };
        let nsam1 = allele_counts1.values().sum::<i64>() as f64;
        let nsam2 = allele_counts2.values().sum::<i64>() as f64;
        for u in unique_allele_labels {
            let c1 = if let Some(c) = allele_counts1.get(&u) {
                *c
            } else {
                0
            };
            let c2 = if let Some(c) = allele_counts2.get(&u) {
                *c
            } else {
                0
            };
            let p1 = (c1 as f64) / nsam1;
            let p2 = (c2 as f64) / nsam2;
            naive += (p1 - p2) * (p1 - p2);
            n += 1;
        }
    }
    naive / (n as f64)
}
