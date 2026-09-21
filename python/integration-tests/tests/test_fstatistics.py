import demes
import msprime
import numpy as np

from hypothesis import given
from hypothesis.strategies import integers

import integration_tests


@given(anc_seed=integers(1, 42000000), mut_seed=integers(1, 42000000))
def test_f2_all_sample_nodes(anc_seed, mut_seed):
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
    counts = integration_tests.counts_from_ts_holder_multi_sample_sets(
        tsholder,
        [
            [i for i in range(ts.num_nodes) if ts.node(i).population == 1],
            [i for i in range(ts.num_nodes) if ts.node(i).population == 2],
        ])
