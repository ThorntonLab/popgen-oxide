import random

import demes
import msprime
import numpy as np
import tskit

from hypothesis import given, reproduce_failure
from hypothesis.strategies import integers

import integration_tests


@given(anc_seed=integers(1, 42000000), mut_seed=integers(1, 42000000))
def test_f2_all_sample_nodes_infinite_sites_mutation(anc_seed, mut_seed):
    yaml = """
    time_units: generations
    demes:
     - name: ancestor
       epochs:
        - start_size: 10000
          end_time: 500
     - name: d0
       ancestors: [ancestor]
       epochs:
        - start_size: 10000
     - name: d1
       ancestors: [ancestor]
       epochs:
        - start_size: 10000
    """
    graph = demes.loads(yaml)
    demog = msprime.Demography.from_demes(graph)
    ts = msprime.sim_ancestry([msprime.SampleSet(100, population=1), msprime.SampleSet(100, population=2)], sequence_length=1000000, demography=demog,
                              recombination_rate=1.2e-8, random_seed=anc_seed)
    ts = msprime.sim_mutations(
        ts, rate=1.2e-8, random_seed=mut_seed, model=msprime.InfiniteSites())
    tsholder = integration_tests.ts_holder_from_tables(ts.tables.copy())
    samples = [
        [i for i in ts.samples(population=1)],
        [i for i in ts.samples(population=2)],
    ]
    counts = integration_tests.counts_from_ts_holder_multi_sample_sets(
        tsholder, samples)
    fstats = integration_tests.fstats(counts)
    f2 = fstats.f2(0, 1)
    f2_py = ts.f2(sample_sets=samples, span_normalise=False)
    for pop in [0, 1]:
        tsdiv = ts.diversity(sample_sets=[samples[pop]], span_normalise=False)
        assert np.isclose(fstats.diversity(pop), tsdiv[0], 1e-10)

    tsdiv = ts.divergence(
        sample_sets=samples, span_normalise=False)
    assert np.isclose(fstats.divergence(0, 1), tsdiv, 1e-10)
    assert np.isclose(f2, f2_py, 1e-10)


@given(anc_seed=integers(1, 42000000), mut_seed=integers(1, 42000000),
       # NOTE: we get failures with smaller sample sizes, which we have to investigate later
       num_sample_nodes_deme1=integers(8, 100),
       num_sample_nodes_deme2=integers(8, 100))
@reproduce_failure('6.155.7', b'AEJc4kEBQQhBCA==')
def test_f2_subset_sample_nodes_infinite_sites_mutation(anc_seed, mut_seed, num_sample_nodes_deme1, num_sample_nodes_deme2):
    yaml = """
    time_units: generations
    demes:
     - name: ancestor
       epochs:
        - start_size: 10000
          end_time: 500
     - name: d0
       ancestors: [ancestor]
       epochs:
        - start_size: 10000
     - name: d1
       ancestors: [ancestor]
       epochs:
        - start_size: 10000
    """
    graph = demes.loads(yaml)
    demog = msprime.Demography.from_demes(graph)
    random.seed(mut_seed)  # ARBITRARY!!
    ts = msprime.sim_ancestry([msprime.SampleSet(100, population=1),
                               msprime.SampleSet(100, population=2)],
                              sequence_length=1000000, demography=demog,
                              recombination_rate=1.2e-8, random_seed=anc_seed)
    ts = msprime.sim_mutations(
        ts, rate=1.2e-8, random_seed=mut_seed, model=msprime.InfiniteSites())
    tsholder = integration_tests.ts_holder_from_tables(ts.tables.copy())
    samples = [
        [i for i in ts.samples(population=1)],
        [i for i in ts.samples(population=2)],
    ]
    random.shuffle(samples[0])
    random.shuffle(samples[1])
    subsamples = [samples[0][:num_sample_nodes_deme1],
                  samples[1][:num_sample_nodes_deme2]]
    counts = integration_tests.counts_from_ts_holder_multi_sample_sets(
        tsholder, subsamples)
    fstats = integration_tests.fstats(counts)
    f2 = fstats.f2(0, 1)
    f2_py = ts.f2(sample_sets=subsamples, span_normalise=False)
    for pop in [0, 1]:
        tsdiv = ts.diversity(
            sample_sets=[subsamples[pop]], span_normalise=False)
        assert np.isclose(fstats.diversity(pop), tsdiv[0], 1e-10)

    tsdiv = ts.divergence(
        sample_sets=subsamples, span_normalise=False)
    assert np.isclose(fstats.divergence(0, 1), tsdiv, 1e-10)

    assert np.isclose(f2, f2_py, 1e-10), f"{f2} {f2_py}"


def test_reciprocal_fixation():
    tables = tskit.TableCollection(10.)
    n0 = tables.nodes.add_row(tskit.NODE_IS_SAMPLE, 0.0)
    n1 = tables.nodes.add_row(tskit.NODE_IS_SAMPLE, 0.0)
    n2 = tables.nodes.add_row(tskit.NODE_IS_SAMPLE, 0.0)
    n3 = tables.nodes.add_row(tskit.NODE_IS_SAMPLE, 0.0)
    _n4 = tables.nodes.add_row(0, 1.)
    n5 = tables.nodes.add_row(0, 1.)
    _n6 = tables.nodes.add_row(0, 2.)

    _ = tables.edges.add_row(0., 10., _n4, n0)
    _ = tables.edges.add_row(0., 10., _n4, n1)
    _ = tables.edges.add_row(0., 10., n5, n2)
    _ = tables.edges.add_row(0., 10., n5, n3)
    _ = tables.edges.add_row(0., 10., _n6, _n4)
    _ = tables.edges.add_row(0., 10., _n6, n5)

    s = tables.sites.add_row(5., ancestral_state='A')
    _ = tables.mutations.add_row(node=n5, site=s, time=1.1, derived_state='G')

    tables.sort()
    tables.build_index()
    ts = tables.tree_sequence()
    tsholder = integration_tests.ts_holder_from_tables(ts.tables.copy())
    samples = [
        [0, 1],
        [2, 3],
    ]
    counts = integration_tests.counts_from_ts_holder_multi_sample_sets(
        tsholder, samples)
    fstats = integration_tests.fstats(counts)
    for pop in [0, 1]:
        tsdiv = ts.diversity(
            sample_sets=[samples[pop]], span_normalise=False)
        assert np.isclose(fstats.diversity(pop), tsdiv[0], 1e-10)

    tsdiv = ts.divergence(
        sample_sets=samples, span_normalise=False)
    assert np.isclose(fstats.divergence(0, 1), tsdiv, 1e-10)
