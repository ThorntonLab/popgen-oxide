import allel
import demes
import msprime
import numpy as np
import pytest

from hypothesis import given
from hypothesis.strategies import integers

import random

import python_integration_tests


def get_allel_thetaw(matrix):
    h = allel.HaplotypeArray(matrix)
    g = h.to_genotypes(ploidy=1)
    ac = g.count_alleles()
    return allel.watterson_theta([i+1 for i in range(h.n_variants)], ac).sum()*(h.n_variants)


@given(anc_seed=integers(1, 42000000), mut_seed=integers(1, 42000000))
def test_thetaw_all_sites(anc_seed, mut_seed):
    yaml = """
    time_units: generations
    demes:
     - name: the_one_deme
       epochs:
        - start_size: 10000
    """
    graph = demes.loads(yaml)
    demog = msprime.Demography.from_demes(graph)
    ts = msprime.sim_ancestry(100, sequence_length=1000000, demography=demog,
                              recombination_rate=1.2e-8, random_seed=anc_seed)
    # NOTE: allele calls a variant site 1 and a monomorphic site 0.
    # We label a site as number of mutations = number of alleles present - 1
    # So, to compare correctness, we need to restrict input to biallelic data
    ts = msprime.sim_mutations(
        ts, rate=1.2e-8, random_seed=mut_seed, discrete_genome=False)

    tsholder = python_integration_tests.ts_holder_from_tables(ts.tables.copy())
    counts = python_integration_tests.counts_from_ts_holder(tsholder)
    thetaw = python_integration_tests.thetaw(counts)
    allel_thetaw = get_allel_thetaw(ts.genotype_matrix())
    assert np.fabs(allel_thetaw-thetaw) <= 1e-9


# @given(anc_seed=integers(1, 42000000), mut_seed=integers(1, 42000000))
# @pytest.mark.xfail(reason="allel does not allow the denominator to vary<F8>l")
# def test_thetaw_all_sites_with_missing_data(anc_seed, mut_seed):
#     yaml = """
#     time_units: generations
#     demes:
#      - name: the_one_deme
#        epochs:
#         - start_size: 10000
#     """
#     graph = demes.loads(yaml)
#     demog = msprime.Demography.from_demes(graph)
#     ts = msprime.sim_ancestry(100, sequence_length=1000000, demography=demog,
#                               recombination_rate=1.2e-8, random_seed=anc_seed)
#     # NOTE: allele calls a variant site 1 and a monomorphic site 0.
#     # We label a site as number of mutations = number of alleles present - 1
#     # So, to compare correctness, we need to restrict input to biallelic data
#     ts = msprime.sim_mutations(
#         ts, rate=1.2e-8, random_seed=mut_seed, discrete_genome=False)
#
#     random.seed(anc_seed)
#     g = ts.genotype_matrix()
#     h = allel.HaplotypeArray(g)
#     samples = [i for i in ts.samples()]
#     for hap in range(h.n_haplotypes):
#         nmissing = random.randint(0, ts.num_samples)
#         assert nmissing <= ts.num_samples
#         mdata = random.sample(samples, nmissing)
#         assert len(h[hap, :]) == ts.num_samples
#         h[hap, :][mdata] = -1
#         if nmissing > 0:
#             assert any([i == -1 for i in h[hap, :]])
#     ac = h.to_genotypes(ploidy=1).count_alleles()
#
#     counts = python_integration_tests.single_sample_set_counts_from_haplotype_matrix(
#         h, ts.num_samples)
#     thetaw = python_integration_tests.thetaw(counts)
#     allel_thetaw = allel.watterson_theta(
#         [i+1 for i in range(h.n_variants)], ac).sum()
#     nvariable = 0
#     for i in ac:
#         if len([j for j in i if j > 0]) > 1:
#             nvariable += 1
#     assert np.fabs(allel_thetaw*nvariable -
#                    thetaw) <= 1e-9, f"{allel_thetaw} * {nvariable} ({h.n_variants}, {ac.n_variants}) !=  {thetaw}"

def test_thetaw_missing_data():
    # Here, the "max" sample size is 5,
    # which is what allele will use for the denominator
    c = [[0, -1, 0, 0, 1],
         [0, 1, 0, 0, 1]]
    h = allel.HaplotypeArray(c)
    assert h.n_variants == 2
    g = h.to_genotypes(ploidy=1)
    ac = g.count_alleles()
    counts = python_integration_tests.single_sample_set_counts_from_haplotype_matrix(
        h, 5)
    thetaw = python_integration_tests.thetaw(counts)
    allel_theta = allel.watterson_theta([1, 2], ac)
    assert np.fabs(allel_theta.sum()*2 - thetaw) >= 1e-9

    # NOTE: the values from allel (for some reason) are "per site"
    assert np.fabs(2 * allel_theta - 2/sum([1/i for i in range(1, 5)])) <= 1e-9
    thetaw_raw = 1/sum([1/i for i in range(1, 4)]) + 1 / \
        sum([1/i for i in range(1, 5)])
    assert np.fabs(thetaw-thetaw_raw) <= 1e-9
