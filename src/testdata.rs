#![allow(dead_code)]

use crate::Count;
use rand::prelude::*;

#[derive(Debug, Clone)]
pub struct GenotypeData {
    data: Box<[Count]>,
}

impl GenotypeData {
    pub fn iter(&self) -> impl Iterator<Item = Count> + '_ {
        self.data.iter().cloned()
    }

    fn make_missing_data(self, rng: &mut StdRng, rate: f64) -> Self {
        assert!((0. ..=1.0).contains(&rate));
        let u = rand_distr::Uniform::new_inclusive(0., 1.0).unwrap();
        let data = self
            .data
            .iter()
            .cloned()
            .map(|c| {
                let mut temp = c;
                for _ in 0..c {
                    if u.sample(rng) < rate {
                        temp -= 1;
                    }
                }
                temp
            })
            .collect::<Vec<_>>()
            .into_boxed_slice();
        Self { data }
    }
}

// A *column* of a VCF file, for example
// Ploidy is implicit
// pub struct Individual {
//     data: Box<[GenotypeData]>,
// }

// A row of a VCF file, for example
// Ploidy is implicit
#[derive(Debug, Clone)]
pub struct Site {
    data: Box<[GenotypeData]>,
}

impl Site {
    fn random(
        rng: &mut StdRng,
        ploidy: usize,
        num_samples: usize,
        allele_frequencies: &[f64],
    ) -> Self {
        let mut data = vec![];
        for _ in 0..num_samples {
            // Collect a multinomial sample
            let temp = multinomial(rng, ploidy, allele_frequencies);
            let temp = temp
                .iter()
                .cloned()
                .map(|i| i64::try_from(i).unwrap())
                .collect::<Vec<_>>()
                .into_boxed_slice();
            let temp = GenotypeData { data: temp };
            data.push(temp)
        }
        assert_eq!(data.len(), num_samples);
        let data = data.into_boxed_slice();
        Self { data }
    }

    fn make_missing(self, rng: &mut StdRng, rate: f64) -> Self {
        let mut x = self;
        x.data = x
            .data
            .iter_mut()
            .map(|x| x.clone().make_missing_data(rng, rate))
            .collect::<Vec<_>>()
            .into_boxed_slice();
        x
    }

    pub fn iter(&self) -> impl Iterator<Item = &GenotypeData> + '_ {
        self.data.iter()
    }
}

#[derive(Clone, Default)]
pub struct RandomSiteOptions {
    missing_data_rate: Option<f64>,
}

// Generate random data at a "site".
// In other words, for a given set of allele frequencies
// at a specific locus, repeatedly sample alleles for a given
// number of samples.
//
// Parameters
//
// * seed: random number seed
// * num_samples: the number of individuals to generate genotypes for
// * ploidy: number of copies of a sites in an individual
// * allele_frequencies: mutations frequencies at the site.
//                       The slice length determines the number of alleles
//                       at the site.
// * options: see RandomSiteOptions.
pub fn random_site(
    seed: u64,
    num_samples: usize,
    ploidy: usize,
    allele_frequencies: &[f64],
    options: Option<RandomSiteOptions>,
) -> Site {
    let options = options.unwrap_or_default();
    let mut rng = StdRng::seed_from_u64(seed);
    let site = Site::random(&mut rng, ploidy, num_samples, allele_frequencies);
    if let Some(rate) = options.missing_data_rate {
        site.make_missing(&mut rng, rate)
    } else {
        site
    }
}

// Helper fns below

// The rand crate has no multinomial function so we will make do with
// rolling our own using the conditional binomial method.
//
// We keep this as a 'vanilla' implementation that doesn't use any
// types defined in the crate.
fn multinomial(rng: &mut StdRng, n: usize, coefficients: &[f64]) -> Box<[usize]> {
    let sum: f64 = coefficients.iter().sum();
    assert!(sum > 0.0);
    let normalized = coefficients
        .iter()
        .cloned()
        .map(|v| v / sum)
        .collect::<Vec<_>>();

    let mut rv = vec![0_usize; normalized.len()];

    let mut nremaining = n;
    let mut cum = 1.0;

    for allele in 0..rv.len() {
        let p = normalized[allele] / cum;
        let b = rand_distr::Binomial::new(nremaining as u64, p).unwrap();
        let ni = b.sample(rng);
        rv[allele] = ni.try_into().unwrap();
        nremaining -= ni as usize;
        cum -= normalized[allele];
        assert!(cum >= 0.0);
        if nremaining == 0 {
            break;
        }
    }
    assert_eq!(nremaining, 0);
    assert_eq!(rv.iter().cloned().sum::<usize>(), n);

    rv.into_boxed_slice()
}

// Tests of our test fns. (A bit meta...)

#[test]
fn test_random_site() {
    let mut rng = StdRng::seed_from_u64(0);
    let s = Site::random(&mut rng, 2, 10, &[0.25, 0.75]);
    s.iter().for_each(|g| {
        assert_eq!(g.iter().sum::<Count>(), 2);
    });

    let c = s.clone().make_missing(&mut rng, 1.0);
    c.iter().for_each(|g| {
        assert!(g.iter().all(|i| i == 0), "{c:?}");
    });

    let c = s.clone().make_missing(&mut rng, 0.0);
    c.iter().for_each(|g| {
        assert_eq!(g.iter().sum::<Count>(), 2);
    });
}
