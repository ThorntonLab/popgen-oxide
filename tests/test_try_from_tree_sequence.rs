// NOTE: these tests require compiling
// with the tskit feature

use varistat::tskit;
use varistat::tskit::prelude::StreamingIterator;

#[cfg(test)]
struct MutationData {
    node: tskit::NodeId,
    time: tskit::Time,
    derived_state: Vec<u8>,
}

#[cfg(test)]
impl MutationData {
    fn new<N, T>(node: N, time: T, derived_state: &str) -> Self
    where
        N: Into<tskit::NodeId>,
        T: Into<tskit::Time>,
    {
        Self {
            node: node.into(),
            time: time.into(),
            derived_state: derived_state.as_bytes().to_owned(),
        }
    }
}

#[cfg(test)]
struct SiteData {
    position: tskit::Position,
    ancestral_state: Vec<u8>,
    mutations: Vec<MutationData>,
}

#[cfg(test)]
impl SiteData {
    fn new<P, M>(position: P, ancestral_state: &str, mutations: M) -> Self
    where
        P: Into<tskit::Position>,
        M: IntoIterator<Item = MutationData>,
    {
        Self {
            position: position.into(),
            ancestral_state: ancestral_state.as_bytes().to_owned(),
            mutations: mutations.into_iter().collect::<Vec<_>>(),
        }
    }
}

#[cfg(test)]
#[derive(Debug)]
struct DerivedCounts {
    count: i64,
    number_of_sites: usize,
}

#[cfg(test)]
#[derive(Debug)]
struct SiteCountContents {
    num_ancestral: i64,
    derived: Vec<DerivedCounts>,
}

#[cfg(test)]
fn validate_site_counts(counts: &varistat::iter::SiteCounts, expected: SiteCountContents) {
    let num_alleles = expected
        .derived
        .iter()
        .map(|s| s.number_of_sites)
        .sum::<usize>()
        + 1;
    let c = counts.counts();
    assert_eq!(c[0], expected.num_ancestral);
    expected.derived.iter().for_each(|d| {
        assert_eq!(
            c.iter().skip(1).filter(|&&i| i == d.count).count(),
            d.number_of_sites
        )
    });
    assert_eq!(c.len(), num_alleles);
}

#[cfg(test)]
mod naive_details {
    use super::SiteCountContents;

    pub fn process_site(
        ts: &tskit::TreeSequence,
        tree: &tskit::Tree,
        current_site: u64,
        current_mutation: u64,
        expected: &mut Vec<SiteCountContents>,
    ) -> u64 {
        let (ancestral_state, mut sample_node_state) =
            set_up_sample_node_states(ts, current_site as i32);
        let current_mutation = update_sample_node_states(
            ts,
            tree,
            current_mutation,
            current_site,
            &mut sample_node_state,
        );
        update_return_value(&ancestral_state, &sample_node_state, expected);
        current_mutation
    }

    struct SampleNodeStates {
        sample_node_state: Vec<Vec<u8>>,
        sample_map: Vec<Option<usize>>,
    }

    impl std::ops::Index<usize> for SampleNodeStates {
        type Output = Vec<u8>;
        fn index(&self, index: usize) -> &Self::Output {
            &self.sample_node_state[self.sample_map[index].unwrap()]
        }
    }

    impl std::ops::IndexMut<usize> for SampleNodeStates {
        fn index_mut(&mut self, index: usize) -> &mut Self::Output {
            &mut self.sample_node_state[self.sample_map[index].unwrap()]
        }
    }

    fn set_up_sample_node_states(
        ts: &tskit::TreeSequence,
        site: i32,
    ) -> (Vec<u8>, SampleNodeStates) {
        let ancestral_state = ts.sites().ancestral_state(site).unwrap().to_owned();
        let sample_node_state =
            vec![ancestral_state.clone(); usize::try_from(ts.num_samples()).unwrap()];
        let mut sample_map = vec![None; ts.nodes().num_rows().try_into().unwrap()];
        for (next_sample, &sample) in ts.sample_nodes().iter().enumerate() {
            sample_map[usize::try_from(sample).unwrap()] = Some(next_sample);
        }
        (
            ancestral_state,
            SampleNodeStates {
                sample_node_state,
                sample_map,
            },
        )
    }

    fn update_sample_node_states(
        ts: &tskit::TreeSequence,
        tree: &tskit::Tree,
        current_mutation: u64,
        current_site: u64,
        sample_node_state: &mut SampleNodeStates,
    ) -> u64 {
        let mut current_mutation = current_mutation;
        while current_mutation < ts.mutations().num_rows()
            && ts.mutations().site(current_mutation as i32).unwrap() == (current_site as i32)
        {
            let node = ts.mutations().node(current_mutation as i32).unwrap();
            for sample in tree.samples(node).unwrap() {
                sample_node_state[usize::try_from(sample).unwrap()] = ts
                    .mutations()
                    .derived_state(current_mutation as i32)
                    .unwrap()
                    .to_owned();
            }
            current_mutation += 1;
        }
        current_mutation
    }

    fn get_unique_alleles(sample_node_state: &SampleNodeStates) -> Vec<Vec<u8>> {
        let mut rv = sample_node_state.sample_node_state.clone();
        rv.sort_unstable();
        rv.dedup();
        rv
    }

    fn count_ancestral_state(
        ancestral_state: &[u8],
        unique_allele_states: &[Vec<u8>],
        sample_node_state: &SampleNodeStates,
    ) -> usize {
        if let Some(a) = unique_allele_states.iter().find(|&u| u == ancestral_state) {
            sample_node_state
                .sample_node_state
                .iter()
                .filter(|&u| u == a)
                .count()
        } else {
            0
        }
    }

    fn count_derived_alleles(
        ancestral_state: &[u8],
        unique_allele_states: &[Vec<u8>],
        sample_node_state: &SampleNodeStates,
    ) -> Vec<usize> {
        // any allele can be present [0, num_samples] times,
        // giving num_samples - 0 + 1 possible values.
        let mut num_derived_counts = vec![0; sample_node_state.sample_node_state.len() + 1];
        unique_allele_states
            .iter()
            .filter(|&u| u != ancestral_state)
            .map(|d| {
                sample_node_state
                    .sample_node_state
                    .iter()
                    .filter(|&u| u == d)
                    .count()
            })
            .for_each(|c| num_derived_counts[c] += 1);
        num_derived_counts
    }

    fn update_return_value(
        ancestral_state: &[u8],
        sample_node_state: &SampleNodeStates,
        expected: &mut Vec<SiteCountContents>,
    ) {
        let unique_allele_states = get_unique_alleles(sample_node_state);
        if unique_allele_states.len() > 1 {
            let mut current_site_count_data = vec![];
            let num_ancestral_allele =
                count_ancestral_state(ancestral_state, &unique_allele_states, sample_node_state);
            let derived_allele_counts =
                count_derived_alleles(ancestral_state, &unique_allele_states, sample_node_state);
            for (i, j) in derived_allele_counts
                .iter()
                .cloned()
                .enumerate()
                .filter(|(_, num_sites)| num_sites > &0)
            {
                current_site_count_data.push(super::DerivedCounts {
                    count: i as i64,
                    number_of_sites: j,
                })
            }
            expected.push(SiteCountContents {
                num_ancestral: num_ancestral_allele as i64,
                derived: current_site_count_data,
            });
        }
    }
}

#[cfg(test)]
fn generate_expected_site_counts_naive(
    ts: &tskit::TreeSequence,
    options: Option<varistat::FromTreeSequenceOptions>,
) -> Vec<SiteCountContents> {
    let mut expected = vec![];
    // We have no code depending on options yet
    assert!(options.is_none());
    let mut current_site = 0_u64;
    let mut current_mutation = 0_u64;
    // We use the EXPENSIVE option of tracking what sample
    // nodes are descendants of each node in each tree.
    let mut tree_iterator = ts.tree_iterator(tskit::TreeFlags::SAMPLE_LISTS).unwrap();
    while let Some(tree) = tree_iterator.next() {
        let (left, right) = tree.interval();
        while current_site < ts.sites().num_rows()
            && ts.sites().position(current_site as i32).unwrap() < left
        {
            current_site += 1;
        }
        while current_site < ts.sites().num_rows()
            && ts.sites().position(current_site as i32).unwrap() < right
        {
            while current_mutation < ts.mutations().num_rows()
                && ts.mutations().site(current_mutation as i32).unwrap() != (current_site as i32)
            {
                current_mutation += 1;
            }
            current_mutation = naive_details::process_site(
                ts,
                tree,
                current_site,
                current_mutation,
                &mut expected,
            );
            current_site += 1;
        }
    }
    expected
}

//   2  <- time is 10
//  ---
//  | |
//  0 1 <- time is 0
#[cfg(test)]
fn make_two_sample_tree() -> tskit::TableCollection {
    let mut tables = tskit::TableCollection::new(100.).unwrap();
    let n0 = tables
        .add_node(
            tskit::NodeFlags::new_sample(),
            0.0,
            tskit::PopulationId::NULL,
            tskit::IndividualId::NULL,
        )
        .unwrap();
    let n1 = tables
        .add_node(
            tskit::NodeFlags::new_sample(),
            0.0,
            tskit::PopulationId::NULL,
            tskit::IndividualId::NULL,
        )
        .unwrap();

    let n2 = tables
        .add_node(
            tskit::NodeFlags::default(),
            10.0,
            tskit::PopulationId::NULL,
            tskit::IndividualId::NULL,
        )
        .unwrap();
    let _ = tables.add_edge(0., 100., n2, n0).unwrap();
    let _ = tables.add_edge(0., 100., n2, n1).unwrap();
    tables
}

//   --6-- <- time 30
//   |   |
//  -5-- | <- time 20
//  |  | |
// -4- | | <- time 10
// | | | |
// 0 1 2 3 <- time 0
#[cfg(test)]
fn make_four_sample_tree() -> tskit::TableCollection {
    let mut tables = tskit::TableCollection::new(100.0).unwrap();

    for _ in 0..4 {
        tables
            .add_node(
                tskit::NodeFlags::new_sample(),
                0.0,
                tskit::PopulationId::NULL,
                tskit::IndividualId::NULL,
            )
            .unwrap();
    }

    for (_, time) in (0..3).zip([10.0, 20.0, 30.0]) {
        tables
            .add_node(
                tskit::NodeFlags::default(),
                time,
                tskit::PopulationId::NULL,
                tskit::IndividualId::NULL,
            )
            .unwrap();
    }

    tables.add_edge(0., 100., 4, 0).unwrap();
    tables.add_edge(0., 100., 4, 1).unwrap();
    tables.add_edge(0., 100., 5, 4).unwrap();
    tables.add_edge(0., 100., 5, 2).unwrap();
    tables.add_edge(0., 100., 6, 5).unwrap();
    tables.add_edge(0., 100., 6, 3).unwrap();

    tables
}

//   --6-- <- time 30
//   |   |
//  -5-- | <- time 20
//  |  | |\
// -4- | | 7 <- time 10
// | | | |
// 0 1 2 3 <- time 0
#[cfg(test)]
fn make_four_sample_tree_with_one_ancient_sample() -> tskit::TableCollection {
    let mut tables = make_four_sample_tree();

    let ancient_sample = tables
        .add_node(tskit::NodeFlags::new_sample(), 10.0, -1, -1)
        .unwrap();
    tables
        .add_edge(0., tables.sequence_length(), 6, ancient_sample)
        .unwrap();

    tables
}

//   --6-- <- time 30
//   |   |
//  -5-- | <- time 20
//  |  | |
// -4- | 7 <- time 10, 7 is a sample
// | | | |
// 0 1 2 3 <- time 0
#[cfg(test)]
fn make_four_sample_tree_with_one_inline_ancient_sample() -> tskit::TableCollection {
    let mut tables = make_four_sample_tree();

    let ancient_sample = tables
        .add_node(tskit::NodeFlags::new_sample(), 10.0, -1, -1)
        .unwrap();
    let mut edges = tskit::EdgeTable::default();

    // copied from make_four_sample_tree
    edges.add_row(0., 100., 4, 0).unwrap();
    edges.add_row(0., 100., 4, 1).unwrap();
    edges.add_row(0., 100., 5, 4).unwrap();
    edges.add_row(0., 100., 5, 2).unwrap();
    edges.add_row(0., 100., 6, 5).unwrap();

    // differences from make_four_sample_tree
    tables
        .add_edge(0., tables.sequence_length(), 6, ancient_sample)
        .unwrap();
    edges
        .add_row(0., tables.sequence_length(), ancient_sample, 3)
        .unwrap();

    tables.set_edges(&edges).unwrap();

    tables
}

//   This tree is repeated 2x.
//   Once on [0, 50)
//   Again on [50, 100)
//   --6-- <- time 30
//   |   |
//  -5-- | <- time 20
//  |  | |
// -4- | | <- time 10
// | | | |
// 0 1 2 3 <- time 0
#[cfg(test)]
fn make_two_identical_four_sample_trees() -> tskit::TableCollection {
    let mut tables = tskit::TableCollection::new(100.0).unwrap();

    for _ in 0..4 {
        tables
            .add_node(
                tskit::NodeFlags::new_sample(),
                0.0,
                tskit::PopulationId::NULL,
                tskit::IndividualId::NULL,
            )
            .unwrap();
    }

    for (_, time) in (0..3).zip([10.0, 20.0, 30.0]) {
        tables
            .add_node(
                tskit::NodeFlags::default(),
                time,
                tskit::PopulationId::NULL,
                tskit::IndividualId::NULL,
            )
            .unwrap();
    }

    tables.add_edge(0., 50., 4, 0).unwrap();
    tables.add_edge(0., 50., 4, 1).unwrap();
    tables.add_edge(0., 50., 5, 4).unwrap();
    tables.add_edge(0., 50., 5, 2).unwrap();
    tables.add_edge(0., 50., 6, 5).unwrap();
    tables.add_edge(0., 50., 6, 3).unwrap();

    tables.add_edge(50., 100., 4, 0).unwrap();
    tables.add_edge(50., 100., 4, 1).unwrap();
    tables.add_edge(50., 100., 5, 4).unwrap();
    tables.add_edge(50., 100., 5, 2).unwrap();
    tables.add_edge(50., 100., 6, 5).unwrap();
    tables.add_edge(50., 100., 6, 3).unwrap();

    tables
}

//   Tree for [0, 50)
//   --6-- <- time 30
//   |   |
//  -5-- | <- time 20
//  |  | |
// -4- | | <- time 10
// | | | |
// 0 1 2 3 <- time 0
// Tree for [50, 100)
//  --6---
//  |  | |
// -5- | |
// | | | |
// 4 | | |
// | | | |
// 0 1 2 3
#[cfg(test)]
fn make_two_different_four_sample_trees() -> tskit::TableCollection {
    let mut tables = tskit::TableCollection::new(100.0).unwrap();

    for _ in 0..4 {
        tables
            .add_node(
                tskit::NodeFlags::new_sample(),
                0.0,
                tskit::PopulationId::NULL,
                tskit::IndividualId::NULL,
            )
            .unwrap();
    }

    for (_, time) in (0..3).zip([10.0, 20.0, 30.0]) {
        tables
            .add_node(
                tskit::NodeFlags::default(),
                time,
                tskit::PopulationId::NULL,
                tskit::IndividualId::NULL,
            )
            .unwrap();
    }

    tables.add_edge(0., 50., 4, 0).unwrap();
    tables.add_edge(0., 50., 4, 1).unwrap();
    tables.add_edge(0., 50., 5, 4).unwrap();
    tables.add_edge(0., 50., 5, 2).unwrap();
    tables.add_edge(0., 50., 6, 5).unwrap();
    tables.add_edge(0., 50., 6, 3).unwrap();

    tables.add_edge(50., 100., 6, 5).unwrap();
    tables.add_edge(50., 100., 6, 2).unwrap();
    tables.add_edge(50., 100., 6, 3).unwrap();
    tables.add_edge(50., 100., 5, 4).unwrap();
    tables.add_edge(50., 100., 5, 1).unwrap();
    tables.add_edge(50., 100., 4, 0).unwrap();

    tables
}

/// 10
/// |\
/// | 9
/// | |\
/// | | 8
/// | | |\
/// | | | 7
/// | | | |\
/// | | | | 6
/// | | | | |\
/// 0 1 2 3 4 5
#[cfg(test)]
fn make_comb_tree() -> tskit::TableCollection {
    let mut tables = tskit::TableCollection::new(10.0).unwrap();
    let n0 = tables
        .add_node(tskit::NodeFlags::new_sample(), 0.0, -1, -1)
        .unwrap();
    let n1 = tables
        .add_node(tskit::NodeFlags::new_sample(), 0.0, -1, -1)
        .unwrap();
    let n2 = tables
        .add_node(tskit::NodeFlags::new_sample(), 0.0, -1, -1)
        .unwrap();
    let n3 = tables
        .add_node(tskit::NodeFlags::new_sample(), 0.0, -1, -1)
        .unwrap();
    let n4 = tables
        .add_node(tskit::NodeFlags::new_sample(), 0.0, -1, -1)
        .unwrap();
    let n5 = tables
        .add_node(tskit::NodeFlags::new_sample(), 0.0, -1, -1)
        .unwrap();
    let n6 = tables
        .add_node(tskit::NodeFlags::default(), 1.0, -1, -1)
        .unwrap();
    let n7 = tables
        .add_node(tskit::NodeFlags::default(), 2.0, -1, -1)
        .unwrap();
    let n8 = tables
        .add_node(tskit::NodeFlags::default(), 3.0, -1, -1)
        .unwrap();
    let n9 = tables
        .add_node(tskit::NodeFlags::default(), 4.0, -1, -1)
        .unwrap();
    let n10 = tables
        .add_node(tskit::NodeFlags::default(), 5.0, -1, -1)
        .unwrap();

    for (p, c) in [
        (n10, n0),
        (n10, n9),
        (n9, n1),
        (n9, n8),
        (n8, n2),
        (n8, n7),
        (n7, n3),
        (n7, n6),
        (n6, n4),
        (n6, n5),
    ] {
        tables.add_edge(0., 10., p, c).unwrap();
    }
    tables
}

#[cfg(test)]
fn add_sites_and_mutations<S>(tables: &mut tskit::TableCollection, data: S)
where
    S: IntoIterator<Item = SiteData>,
{
    for s in data {
        let site = tables
            .add_site(s.position, Some(&s.ancestral_state))
            .unwrap();
        for m in &s.mutations {
            tables
                .add_mutation(
                    site,
                    m.node,
                    tskit::MutationId::NULL,
                    m.time,
                    Some(&m.derived_state),
                )
                .unwrap();
        }
    }
}

#[cfg(test)]
fn make_test_data<F, S>(make_tables: F, data: S) -> tskit::TreeSequence
where
    F: Fn() -> tskit::TableCollection,
    S: IntoIterator<Item = SiteData>,
{
    let mut tables = make_tables();
    add_sites_and_mutations(&mut tables, data);
    tables
        .full_sort(tskit::TableSortOptions::default())
        .unwrap();
    tables.build_index().unwrap();
    assert!(!tables.as_mut_ptr().is_null());
    tables
        .compute_mutation_parents(tskit::MutationParentsFlags::default())
        .unwrap();
    tables
        .tree_sequence(tskit::TreeSequenceFlags::default())
        .unwrap()
}

#[cfg(test)]
fn generate_counts_and_validate(
    ts: &tskit::TreeSequence,
    options: Option<varistat::FromTreeSequenceOptions>,
) {
    let counts = varistat::MultiSiteCounts::try_from_tree_sequence(ts, None).unwrap();
    let expected = generate_expected_site_counts_naive(ts, options);
    assert_eq!(counts.len(), expected.len());
    for (obs, exp) in counts.iter().zip(expected.into_iter()) {
        validate_site_counts(&obs, exp);
    }
}

#[test]
fn test_0() {
    let ts = make_test_data(make_two_sample_tree, vec![]);
    generate_counts_and_validate(&ts, None);
}

#[test]
fn test_1() {
    let ts = make_test_data(
        make_two_sample_tree,
        vec![SiteData::new(
            5.0,
            "A",
            vec![MutationData::new(1, 1.0, "G")],
        )],
    );
    generate_counts_and_validate(&ts, None);
}

#[test]
fn test_2() {
    let ts = make_test_data(
        make_two_sample_tree,
        vec![
            SiteData::new(5.0, "A", vec![MutationData::new(1, 1.0, "G")]),
            SiteData::new(
                6.0,
                "T",
                vec![
                    MutationData::new(1, 1.0, "G"),
                    MutationData::new(0, 1.0, "C"),
                ],
            ),
        ],
    );
    generate_counts_and_validate(&ts, None);
}

#[test]
fn test_3() {
    let ts = make_test_data(
        make_four_sample_tree,
        vec![SiteData::new(
            5.,
            "T",
            vec![
                MutationData::new(5, 21.0, "G"),
                MutationData::new(2, 10., "T"),
            ],
        )],
    );
    generate_counts_and_validate(&ts, None);
}

#[test]
fn test_4() {
    let ts = make_test_data(
        make_four_sample_tree,
        vec![SiteData::new(
            5.,
            "G",
            vec![
                MutationData::new(5, 21.0, "A"),
                MutationData::new(4, 10.1, "G"),
                MutationData::new(1, 0.1, "C"),
            ],
        )],
    );
    generate_counts_and_validate(&ts, None);
}

#[test]
fn test_5() {
    let site0 = SiteData::new(
        60.0,
        "G",
        vec![
            MutationData::new(5, 20.1, "A"),
            MutationData::new(4, 10.1, "G"),
            MutationData::new(1, 0.1, "C"),
        ],
    );
    let site1 = SiteData::new(
        40.0,
        "G",
        vec![
            MutationData::new(5, 20.1, "A"),
            MutationData::new(4, 10.1, "G"),
            MutationData::new(1, 0.1, "C"),
        ],
    );
    let ts = make_test_data(make_two_identical_four_sample_trees, vec![site0, site1]);
    generate_counts_and_validate(&ts, None);
}

#[test]
fn test_6() {
    let site0 = SiteData::new(
        60.0,
        "G",
        vec![
            MutationData::new(5, 20.1, "A"),
            MutationData::new(4, 10.1, "G"),
            MutationData::new(1, 0.1, "C"),
        ],
    );
    let site1 = SiteData::new(
        40.0,
        "T",
        vec![
            MutationData::new(3, 20.1, "T"),
            MutationData::new(2, 0.1, "G"),
            MutationData::new(4, 10.1, "G"),
        ],
    );
    let ts = make_test_data(make_two_identical_four_sample_trees, vec![site0, site1]);
    generate_counts_and_validate(&ts, None);
}

#[test]
fn test_7() {
    let site0 = SiteData::new(
        60.0,
        "G",
        vec![
            MutationData::new(5, 20.1, "A"),
            MutationData::new(4, 10.1, "G"),
            MutationData::new(1, 0.1, "C"),
        ],
    );
    let site1 = SiteData::new(
        40.0,
        "T",
        vec![
            MutationData::new(3, 20.1, "T"),
            MutationData::new(2, 0.1, "G"),
            MutationData::new(4, 10.1, "G"),
        ],
    );
    let ts = make_test_data(make_two_different_four_sample_trees, vec![site0, site1]);
    generate_counts_and_validate(&ts, None);
}

#[test]
fn test_8() {
    let site0 = SiteData::new(60.0, "G", vec![]);
    let ts = make_test_data(make_two_different_four_sample_trees, vec![site0]);
    generate_counts_and_validate(&ts, None);
}

#[test]
fn test_9() {
    let site0 = SiteData::new(
        60.0,
        "G",
        vec![
            MutationData::new(5, 20.1, "A"),
            MutationData::new(3, 19.1, "A"),
        ],
    );
    let ts = make_test_data(make_two_different_four_sample_trees, vec![site0]);
    generate_counts_and_validate(&ts, None);
}

#[test]
#[should_panic]
fn test_10_anc_state_missing() {
    let site0 = SiteData::new(
        60.0,
        "",
        vec![
            MutationData::new(5, 20.1, "A"),
            MutationData::new(3, 19.1, "A"),
        ],
    );
    let ts = make_test_data(make_two_different_four_sample_trees, vec![site0]);
    generate_counts_and_validate(&ts, None);
}

#[test]
#[should_panic]
fn test_10_der_state_missing() {
    let site0 = SiteData::new(
        60.0,
        "G",
        vec![
            MutationData::new(5, 20.1, ""),
            MutationData::new(3, 19.1, "A"),
        ],
    );
    let ts = make_test_data(make_two_different_four_sample_trees, vec![site0]);
    generate_counts_and_validate(&ts, None);
}

#[test]
fn test_11() {
    let ts = make_test_data(
        make_four_sample_tree,
        vec![SiteData::new(
            5.,
            "T",
            vec![
                MutationData::new(5, 21.0, "A"),
                MutationData::new(4, 10.1, "G"),
            ],
        )],
    );
    generate_counts_and_validate(&ts, None);
}

#[test]
fn test_12() {
    let ts = make_test_data(
        make_comb_tree,
        vec![SiteData::new(
            5.0,
            "0",
            vec![
                MutationData::new(8, 3.1, "1"),
                MutationData::new(7, 2.1, "2"),
                MutationData::new(6, 1.1, "3"),
                MutationData::new(5, 0.1, "4"),
                MutationData::new(9, 4.1, "5"),
            ],
        )],
    );
    generate_counts_and_validate(&ts, None);
}

#[cfg(test)]
mod with_ancient_samples {
    use super::*;

    #[test]
    fn test0() {
        let ts = super::make_test_data(
            make_four_sample_tree_with_one_ancient_sample,
            vec![SiteData::new(
                5.,
                "G",
                vec![
                    MutationData::new(5, 21.0, "A"),
                    MutationData::new(4, 10.1, "G"),
                    MutationData::new(1, 0.1, "C"),
                ],
            )],
        );
        generate_counts_and_validate(&ts, None);
    }

    #[test]
    fn test1() {
        let ts = super::make_test_data(
            make_four_sample_tree_with_one_ancient_sample,
            vec![SiteData::new(
                5.,
                "G",
                vec![
                    MutationData::new(5, 21.0, "A"),
                    MutationData::new(4, 10.1, "G"),
                    MutationData::new(7, 10.1, "C"),
                    MutationData::new(1, 0.1, "C"),
                ],
            )],
        );
        generate_counts_and_validate(&ts, None);
    }

    #[test]
    fn test2() {
        let ts = super::make_test_data(
            make_four_sample_tree_with_one_ancient_sample,
            vec![SiteData::new(
                5.,
                "G",
                vec![
                    MutationData::new(5, 21.0, "A"),
                    MutationData::new(4, 10.1, "G"),
                    MutationData::new(7, 10.1, "A"),
                    MutationData::new(1, 0.1, "C"),
                ],
            )],
        );
        generate_counts_and_validate(&ts, None);
    }

    #[test]
    fn test3() {
        let ts = super::make_test_data(
            make_four_sample_tree_with_one_inline_ancient_sample,
            vec![SiteData::new(
                5.,
                "G",
                vec![
                    MutationData::new(5, 21.0, "A"),
                    MutationData::new(4, 10.1, "G"),
                    MutationData::new(7, 10.1, "A"),
                    MutationData::new(1, 0.1, "C"),
                ],
            )],
        );
        generate_counts_and_validate(&ts, None);
    }
}
