import demes
import msprime
import numpy as np
import pytest

from hypothesis import given
from hypothesis.strategies import integers

import integration_tests


@given(anc_seed=integers(1, 42000000), mut_seed=integers(1, 42000000))
def test_diversity_all_sites(anc_seed, mut_seed):
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
    ts = msprime.sim_mutations(ts, rate=1.2e-8, random_seed=mut_seed)
    tsdiv = ts.diversity(
        sample_sets=[i for i in ts.samples()], span_normalise=False)

    tsholder = integration_tests.ts_holder_from_tables(ts.tables.copy())
    counts = integration_tests.counts_from_ts_holder(tsholder)
    div = integration_tests.diversity(counts)
    assert np.fabs(tsdiv-div) <= 1e-9


@given(anc_seed=integers(1, 42000000), mut_seed=integers(1, 42000000), np_seed=integers(0, 42000000))
def test_diversity_all_sites_subsample_nodes(anc_seed, mut_seed, np_seed):
    yaml = """
    time_units: generations
    demes:
     - name: the_one_deme
       epochs:
        - start_size: 10000
    """
    np.random.seed(np_seed)
    graph = demes.loads(yaml)
    demog = msprime.Demography.from_demes(graph)
    ts = msprime.sim_ancestry(100, sequence_length=1000000, demography=demog,
                              recombination_rate=1.2e-8, random_seed=anc_seed)
    ts = msprime.sim_mutations(ts, rate=1.2e-8, random_seed=mut_seed)

    all_samples = [i for i in ts.samples()]
    samples = np.random.choice(all_samples, 50, replace=False)
    tsdiv = ts.diversity(
        sample_sets=samples, span_normalise=False)
    tsholder = integration_tests.ts_holder_from_tables(ts.tables.copy())
    counts = integration_tests.counts_from_ts_holder_single_sample_set(
        tsholder, samples)
    div = integration_tests.diversity(counts)
    assert np.fabs(tsdiv-div) <= 1e-9


@given(anc_seed=integers(1, 42000000), mut_seed=integers(1, 42000000), np_seed=integers(0, 42000000))
def test_diversity_all_sites_one_ancient_sample_set(anc_seed, mut_seed, np_seed):
    yaml = """
    time_units: generations
    demes:
     - name: the_one_deme
       epochs:
        - start_size: 10000
    """
    np.random.seed(np_seed)
    graph = demes.loads(yaml)
    demog = msprime.Demography.from_demes(graph)
    ts = msprime.sim_ancestry(samples=[msprime.SampleSet(100), msprime.SampleSet(53, time=214)], sequence_length=1000000, demography=demog,
                              recombination_rate=1.2e-8, random_seed=anc_seed)
    ts = msprime.sim_mutations(ts, rate=1.2e-8, random_seed=mut_seed)

    samples = [i for i in ts.samples()]
    tsdiv = ts.diversity(
        sample_sets=samples, span_normalise=False)
    tsholder = integration_tests.ts_holder_from_tables(ts.tables.copy())
    counts = integration_tests.counts_from_ts_holder_single_sample_set(
        tsholder, samples)
    div = integration_tests.diversity(counts)
    assert np.fabs(tsdiv-div) <= 1e-9


@given(anc_seed=integers(1, 42000000), mut_seed=integers(1, 42000000), np_seed=integers(0, 42000000), num_windows=integers(2, 15))
def test_diversity_all_sites_windows(anc_seed, mut_seed, np_seed, num_windows):
    yaml = """
    time_units: generations
    demes:
     - name: the_one_deme
       epochs:
        - start_size: 10000
    """
    np.random.seed(np_seed)
    graph = demes.loads(yaml)
    demog = msprime.Demography.from_demes(graph)
    ts = msprime.sim_ancestry(samples=[msprime.SampleSet(100)], sequence_length=1000000, demography=demog,
                              recombination_rate=1.2e-8, random_seed=anc_seed)
    ts = msprime.sim_mutations(ts, rate=1.2e-8, random_seed=mut_seed)

    samples = [i for i in ts.samples()]
    windows = np.linspace(0.0, ts.sequence_length, num_windows, endpoint=True)
    tsdiv = ts.diversity(
        sample_sets=samples, span_normalise=False, windows=windows)
    tsholder = integration_tests.ts_holder_from_tables(ts.tables.copy())
    windowed_counts = integration_tests.counts_from_ts_holder_windowed(
        tsholder, samples, windows)
    div = windowed_counts.diversity()
    assert len(tsdiv) == len(div)
    for i, j in zip(tsdiv, div):
        assert np.fabs(i-j) <= 1e-9
