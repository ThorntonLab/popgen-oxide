import demes
import msprime
import numpy as np

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


# This test hits our API harder by not requiring that windows cover the entire
# length of the "genome" by skipping the first and last window
# TODO: write a custom hypothesis strategy for "number of subwindows" and
# we an use numpy to subset the input windows...
@given(anc_seed=integers(1, 42000000), mut_seed=integers(1, 42000000), np_seed=integers(0, 42000000), num_windows=integers(4, 15))
def test_diversity_windows(anc_seed, mut_seed, np_seed, num_windows):
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
    window_subset = windows[1:len(windows)-1]
    windows_tupled = [(i, j)
                      for i, j in zip(window_subset[0:], window_subset[1:])]
    assert len(window_subset) > 0, f"{window_subset}"
    windowed_counts = integration_tests.counts_from_ts_holder_windowed_general(
        tsholder, samples, windows_tupled)
    div = windowed_counts.diversity()
    assert len(tsdiv) - 2 == len(div)
    for i, j in zip(tsdiv[1:len(tsdiv)-1], div):
        assert np.fabs(i-j) <= 1e-9

    if len(window_subset) > 2:
        windowed_counts = integration_tests.counts_from_ts_holder_windowed_general(
            tsholder, samples, windows_tupled[::2])
        div = windowed_counts.diversity()
        for i, j in zip(tsdiv[1:len(tsdiv)-1][::2], div):
            assert np.fabs(i-j) <= 1e-9

    if len(windows) > 3:
        windows_tupled = [(i, j)
                          for i, j in zip(windows[0:], windows[1:])]
        assert len(windows_tupled) == len(windows) - 1
        windowed_counts = integration_tests.counts_from_ts_holder_windowed_general(
            tsholder, samples, windows_tupled[::2])
        div = windowed_counts.diversity()
        for w, d in zip(windows_tupled[::2], div):
            found = False
            for i, j in zip(windows_tupled, tsdiv):
                if i == w:
                    found = True
                    assert np.fabs(
                        d-j) <= 1e-9, f"{w} {i} {d} != {j}, {tsdiv} -> {tsdiv[::2]} vs {div}"
                    break
            assert found is True
