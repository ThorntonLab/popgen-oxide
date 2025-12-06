#![allow(dead_code)]

use crate::MultiSiteCounts;
use rand::prelude::*;

#[derive(Debug, Clone)]
pub struct GenotypeData {
    data: Box<[Option<usize>]>,
}

impl GenotypeData {
    pub fn iter(&self) -> impl Iterator<Item = Option<usize>> + '_ {
        self.data.iter().cloned()
    }

    fn make_missing_data(self, rng: &mut StdRng, rate: f64) -> Self {
        assert!((0. ..=1.0).contains(&rate));
        let u = rand_distr::Uniform::new_inclusive(0., 1.0).unwrap();
        let data = self
            .data
            .iter()
            .cloned()
            .map(|c| if u.sample(rng) < rate { None } else { c })
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
            let mut genotype = vec![];
            for (i, &t) in temp.iter().enumerate() {
                for _ in 0..t {
                    genotype.push(Some(i))
                }
            }
            for _ in genotype.len()..ploidy {
                genotype.push(None);
            }
            assert_eq!(genotype.len(), ploidy);
            let temp = GenotypeData {
                data: genotype.into_boxed_slice(),
            };
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
            .iter()
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
    pub(super) missing_data_rate: Option<f64>,
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
    let mut rng = StdRng::seed_from_u64(seed);
    random_site_rng(num_samples, ploidy, allele_frequencies, options, &mut rng)
}

pub fn random_site_rng(
    num_samples: usize,
    ploidy: usize,
    allele_frequencies: &[f64],
    options: Option<RandomSiteOptions>,
    rng: &mut StdRng,
) -> Site {
    let options = options.unwrap_or_default();
    let site = Site::random(rng, ploidy, num_samples, allele_frequencies);
    if let Some(rate) = options.missing_data_rate {
        site.make_missing(rng, rate)
    } else {
        site
    }
}

pub fn single_pop_counts<'s>(sites: &'s mut dyn Iterator<Item = &'s Site>) -> MultiSiteCounts {
    let mut mcounts = MultiSiteCounts::default();
    for s in sites {
        // The number of INDIVIDUAL genotypes at this site
        let num_samples = s.iter().count();
        if let Some(max_allele_id) = s
            .iter()
            .flat_map(|i| i.iter())
            .flatten()
            .max_by(|i, j| i.cmp(j))
        {
            let ploidy = s.iter().take(1).flat_map(|i| i.iter()).count();
            assert!(ploidy > 0);
            let mut counts = vec![0; max_allele_id + 1];
            s.iter()
                .flat_map(|g| g.iter())
                .flatten()
                .for_each(|i| counts[i] += 1);
            // Probably overly strict but we might as well
            // never let anything invalid slide through test functions
            let ploidy = i32::try_from(ploidy).unwrap();
            let num_samples = i32::try_from(num_samples).unwrap();
            let total_alleles = ploidy.checked_mul(num_samples).unwrap();
            mcounts
                .add_site_from_counts(&counts, total_alleles)
                .unwrap();
        }
    }
    mcounts
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
    assert_eq!(rv.iter().sum::<usize>(), n);

    rv.into_boxed_slice()
}

/// A factory which yields, via [`.iter()`](Self::iter), iterator(s) yielding `Vec<f64>`s which:
/// * Each have length `num_alleles`.
/// * Each yielded `Vec` will have a different index uniquely set to `1.0`.
///   All other elements will be `0.0`.
///
/// The elements yielded from one iterator represent every possible allele frequency distribution for a single *fixed mutation* given that there are `num_alleles` total alleles.
#[derive(Copy, Clone, Debug)]
pub struct FixedMutationIterator {
    num_alleles: usize,
}

impl FixedMutationIterator {
    pub fn new(num_alleles: usize) -> Self {
        assert!(num_alleles > 0);
        Self { num_alleles }
    }

    pub fn iter(&self) -> Box<dyn Iterator<Item = Vec<f64>>> {
        Box::new(FixedMutationIteratorDetail {
            num_alleles: self.num_alleles,
            current_allele: 0,
        })
    }
}

struct FixedMutationIteratorDetail {
    num_alleles: usize,
    current_allele: usize,
}

impl Iterator for FixedMutationIteratorDetail {
    type Item = Vec<f64>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.current_allele < self.num_alleles {
            let mut rv = vec![0.0; self.num_alleles];
            rv[self.current_allele] = 1.0;
            self.current_allele += 1;
            Some(rv)
        } else {
            None
        }
    }
}

pub struct RandomAlleleFreqIterator {
    seed: u64,
    coefficients: Vec<f64>,
}

impl RandomAlleleFreqIterator {
    pub fn new(seed: u64, coefficients: &[f64]) -> Self {
        // We will generally just unwrap Err from rand/rand_distr,
        // but this is one we can catch quite early.
        assert!(coefficients.len() > 1);
        Self {
            seed,
            coefficients: coefficients.to_owned(),
        }
    }

    /// WARNING: return value never returns None, so
    /// caller must use \.take(usize) in order
    /// for things to end.
    pub fn iter(&self) -> Box<dyn Iterator<Item = Vec<f64>>> {
        Box::new(RandomAlleleFreqIteratorDetail::new(
            self.seed,
            &self.coefficients,
        ))
    }
}

// This type is effectively rolling its
// own implementation of rand_distr::Dirichlet.
// We need this b/c the rand_distr version
// is const generic but we require variable
// length outputs at run time.
//
// We do not bother with beta densities
// and just accept the slight efficiency
// loss in some cases. (See the rand_distr
// impl of sample for details.)
struct RandomAlleleFreqIteratorDetail {
    rng: StdRng,
    gammas: Vec<rand_distr::Gamma<f64>>,
}

impl RandomAlleleFreqIteratorDetail {
    fn new(seed: u64, coefficients: &[f64]) -> Self {
        let rng = rand::rngs::StdRng::seed_from_u64(seed);
        let gammas = coefficients
            .iter()
            .cloned()
            .map(|f| rand_distr::Gamma::new(f, 1.0).unwrap())
            .collect::<Vec<_>>();
        Self { rng, gammas }
    }

    // Copy of the rand_distr impl of sample
    fn sample(&mut self) -> Option<Vec<f64>> {
        let mut rv = vec![0.0; self.gammas.len()];
        let mut sum = 0.0;
        for (i, g) in rv.iter_mut().zip(self.gammas.iter()) {
            *i = g.sample(&mut self.rng);
            sum += *i;
        }
        let x = 1. / sum;
        for i in rv.iter_mut() {
            *i *= x;
        }
        Some(rv)
    }
}

impl Iterator for RandomAlleleFreqIteratorDetail {
    type Item = Vec<f64>;

    fn next(&mut self) -> Option<Self::Item> {
        self.sample()
    }
}

// Tests of our test fns. (A bit meta...)

#[test]
fn test_random_site() {
    let mut rng = StdRng::seed_from_u64(0);
    let s = Site::random(&mut rng, 2, 10, &[0.25, 0.75]);
    s.iter().for_each(|g| {
        assert_eq!(g.iter().count(), 2);
    });

    let c = s.clone().make_missing(&mut rng, 1.0);
    c.iter().for_each(|g| {
        assert!(g.iter().all(|i| i.is_none()), "{c:?}");
    });

    let c = s.clone().make_missing(&mut rng, 0.0);
    c.iter().for_each(|g| {
        assert_eq!(g.iter().filter(|i| i.is_some()).count(), 2);
    });
}

#[test]
fn test_fixed_mutation_iterator() {
    let f = FixedMutationIterator::new(4);
    assert_eq!(f.iter().count(), 4);
    for (i, v) in f.iter().enumerate() {
        assert_eq!(v[i], 1.);
        assert_eq!(v.iter().cloned().filter(|&vi| vi == 1.).count(), 1);
        assert_eq!(
            v.iter().cloned().filter(|&vi| vi == 0.).count(),
            v.len() - 1
        );
    }
}

#[test]
fn test_rand_allele_freq_iter() {
    for coefficients in [vec![1., 1., 1.], vec![f64::EPSILON, 1.]] {
        let f = RandomAlleleFreqIterator::new(101, &coefficients);
        for sample in f.iter().take(100) {
            assert_eq!(sample.len(), coefficients.len());
            assert!((sample.iter().sum::<f64>() - 1.).abs() <= 1e-8);
        }
    }
}
