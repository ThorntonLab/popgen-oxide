use crate::testdata::GenotypeData;
use crate::testdata::Site;

pub fn pi_site(genotypes: &mut dyn Iterator<Item = GenotypeData>) -> f64 {
    let mut num_differences = 0_i64;
    let mut num_comparisons = 0_i64;
    let g = genotypes.collect::<Vec<_>>();
    let num_alleles = g[0].iter().count();
    let mut temp = vec![0; num_alleles];
    for (i, j) in g.iter().enumerate() {
        for (x, y) in j.iter().enumerate() {
            temp[x] += y;
        }
        for k in g.iter().skip(i) {
            for b in j.iter() {
                for c in k.iter() {
                    num_differences += (b - c).abs();
                    num_comparisons += b + c;
                }
            }
        }
    }
    let sum: i64 = temp.iter().sum();
    let mut ttl = 0.0;
    for &i in &temp {
        let x = i as f64;
        let n = sum as f64;
        ttl += (x / n) * ((x - 1.) / (n - 1.));
    }
    ttl = 1. - ttl;
    println!(
        "{temp:?} {ttl} {num_differences}, {num_comparisons} {}",
        (num_differences as f64 / num_comparisons as f64)
    );
    num_differences as f64 / num_comparisons as f64
}

pub fn pi<'s>(sites: &'s mut dyn Iterator<Item = &'s mut Site>) -> f64 {
    sites.map(|s| pi_site(&mut s.iter().cloned())).sum::<f64>()
}
