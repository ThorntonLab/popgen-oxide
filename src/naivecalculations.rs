use crate::testdata::GenotypeData;
use crate::testdata::Site;

pub fn pi_site(genotypes: &mut dyn Iterator<Item = GenotypeData>) -> f64 {
    let mut num_differences = 0_i64;
    let mut num_comparisons = 0_i64;
    let g = genotypes.collect::<Vec<_>>();
    for (i, j) in g.iter().enumerate() {
        for k in g.iter().skip(i) {
            for b in j.iter() {
                for c in k.iter() {
                    num_differences += (b - c).abs();
                    num_comparisons += 1;
                }
            }
        }
    }
    println!(
        "{num_differences}, {num_comparisons} {}",
        (num_differences as f64 / num_comparisons as f64)
    );
    num_differences as f64 / num_comparisons as f64
}

pub fn pi<'s>(sites: &'s mut dyn Iterator<Item = &'s mut Site>) -> f64 {
    sites.map(|s| pi_site(&mut s.iter().cloned())).sum::<f64>()
}
