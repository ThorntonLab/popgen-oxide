use crate::testdata::GenotypeData;
use crate::testdata::Site;

fn pi_site(genotypes: &mut dyn Iterator<Item = GenotypeData>) -> f64 {
    let mut num_differences = 0_i64;
    let mut num_comparisons = 0_i64;
    let mut temp = vec![];
    for gi in genotypes {
        for (i, j) in gi.iter().enumerate() {
            for _ in 0..j {
                temp.push(i);
            }
        }
    }
    for (i, j) in temp.iter().enumerate() {
        for k in temp.iter().skip(i + 1) {
            if j != k {
                num_differences += 1;
            }
            num_comparisons += 1;
        }
    }
    assert_eq!(num_comparisons as usize, temp.len() * (temp.len() - 1) / 2);
    num_differences as f64 / num_comparisons as f64
}

// O(N^2) implementation of the Nei/Tajima diversity measure.
pub fn pi<'s>(sites: &'s mut dyn Iterator<Item = &'s mut Site>) -> f64 {
    sites.map(|s| pi_site(&mut s.iter().cloned())).sum::<f64>()
}
