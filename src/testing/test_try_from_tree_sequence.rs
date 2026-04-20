// NOTE: these tests require compiling
// with the tskit feature

use tskit::prelude::StreamingIterator;

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
fn validate_site_counts(counts: &crate::iter::SiteCounts, expected: SiteCountContents) {
    let num_alleles = expected
        .derived
        .iter()
        .map(|s| s.number_of_sites)
        .sum::<usize>()
        + 1;
    let c = counts.counts();
    // FIXME: the following assertion is wrong
    // when playing w/a subset of nodes.
    assert_eq!(c[0], expected.num_ancestral, "{counts:?} {expected:?}");
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

    struct MutationInfo {
        node: tskit::NodeId,
        derived_state: Vec<u8>,
    }

    pub fn process_site(
        ts: &tskit::TreeSequence,
        tree: &tskit::Tree,
        current_site: u64,
        current_mutation: u64,
        focal_nodes: &[bool],
        expected: &mut Vec<SiteCountContents>,
    ) -> u64 {
        let ancestral_state = ts
            .sites()
            .ancestral_state(current_site as i32)
            .unwrap()
            .to_owned();
        // mutations_at_site will contain each mutation at this site,
        // ordered from present to past.
        let (mutations_at_site, current_mutation) = {
            let mut current_mutation = current_mutation;
            let mut mutations_at_site = vec![];
            while current_mutation < ts.mutations().num_rows()
                && ts.mutations().site(current_mutation as i32).unwrap() == (current_site as i32)
            {
                let node = ts.mutations().node(current_mutation as i32).unwrap();
                let derived_state = ts
                    .mutations()
                    .derived_state(current_mutation as i32)
                    .unwrap()
                    .to_owned();
                mutations_at_site.push(MutationInfo {
                    node,
                    derived_state,
                });
                current_mutation += 1;
            }
            mutations_at_site.reverse();
            (mutations_at_site, current_mutation)
        };

        // Use brute-force searching of mutations_at_site to set
        // the state of each node.
        let node_state = {
            let mut node_state = vec![ancestral_state.clone(); ts.nodes().num_rows().as_usize()];
            let parent = tree.parent_array();
            for node in focal_nodes
                .iter()
                .cloned()
                .enumerate()
                .filter_map(|(n, f)| if f { Some(n) } else { None })
            {
                let mut p = tskit::NodeId::from(node as i32);
                while p != tskit::NodeId::NULL {
                    if let Some(index) = mutations_at_site.iter().position(|t| t.node == p) {
                        node_state[node] = mutations_at_site[index].derived_state.clone();
                        break;
                    }
                    p = parent[p.as_usize()];
                }
            }
            node_state
        };

        // extract focal node state
        let focal_node_state = focal_nodes
            .iter()
            .cloned()
            .enumerate()
            .filter_map(|(n, f)| if f { Some(node_state[n].clone()) } else { None })
            .collect::<Vec<_>>();
        let num_ancestral = focal_node_state
            .iter()
            .filter(|i| i == &&ancestral_state)
            .count() as i64;
        let unique_focal_node_derived_states = {
            let mut unique_focal_node_derived_states = focal_node_state
                .iter()
                .filter(|i| i != &&ancestral_state)
                .cloned()
                .collect::<Vec<_>>();
            unique_focal_node_derived_states.sort_unstable();
            unique_focal_node_derived_states.dedup();
            unique_focal_node_derived_states
        };
        let num_derived_counts = {
            let mut num_derived_counts = vec![0; focal_node_state.len() + 1];
            for u in unique_focal_node_derived_states {
                let count = focal_node_state.iter().filter(|&i| i == &u).count();
                num_derived_counts[count] += 1;
            }
            num_derived_counts
        };
        let derived = num_derived_counts
            .iter()
            .cloned()
            .enumerate()
            .filter(|(_, num_sites)| num_sites > &0)
            .map(|(count, number_of_sites)| super::DerivedCounts {
                count: count as i64,
                number_of_sites,
            })
            .collect::<Vec<_>>();
        if !derived.is_empty() {
            expected.push(SiteCountContents {
                num_ancestral,
                derived,
            });
        }
        current_mutation
    }
}

#[cfg(test)]
fn generate_expected_site_counts_naive(
    ts: &tskit::TreeSequence,
    focal_nodes: &[bool],
    _options: Option<&popgen::FromTreeSequenceOptions>,
) -> Vec<SiteCountContents> {
    let mut expected = vec![];
    let mut current_site = 0_u64;
    let mut current_mutation = 0_u64;
    let mut tree_iterator = ts.tree_iterator(0).unwrap();
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
                focal_nodes,
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
    options: Option<&popgen::FromTreeSequenceOptions>,
) {
    let counts = popgen::MultiSiteCounts::try_from_tree_sequence(ts, options).unwrap();
    println!("{counts:?}");
    // Instead of relying on the internal node sample-ness status,
    // we define our set of "sample/focal" nodes externally from
    // the tree sequence.
    let focal_nodes = {
        let mut focal_nodes = vec![false; ts.nodes().num_rows().as_usize()];
        if let Some(opts) = &options {
            if let Some(samples) = &opts.samples {
                match samples {
                    popgen::TskitSamplesList::Node(nodes) => {
                        for n in nodes.iter().cloned().map(|i| i.as_usize()) {
                            assert!(!focal_nodes[n]);
                            focal_nodes[n] = true;
                        }
                    }
                    popgen::TskitSamplesList::Individual(_) => todo!(),
                }
            }
        } else {
            // Get the nodes from the tree sequence sample map
            for n in ts.nodes_iter().filter_map(|n| {
                if n.flags.is_sample() {
                    Some(n.id)
                } else {
                    None
                }
            }) {
                assert!(!focal_nodes[n.as_usize()]);
                focal_nodes[n.as_usize()] = true;
            }
        }
        focal_nodes
    };
    let expected = generate_expected_site_counts_naive(ts, &focal_nodes, options);
    assert_eq!(counts.len(), expected.len(), "{counts:?}, {expected:?}");
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

#[test]
fn test_13() {
    for (mut_node, mut_time) in [(5, 0.1), (6, 1.1), (7, 2.1), (8, 3.1), (9, 4.1)] {
        let ts = make_test_data(
            make_comb_tree,
            vec![SiteData::new(
                5.0,
                "0",
                vec![MutationData::new(mut_node, mut_time, "1")],
            )],
        );
        generate_counts_and_validate(&ts, None);
        let samples = ts.sample_nodes();
        for x in 2..samples.len() - 1 {
            let options = popgen::FromTreeSequenceOptions {
                samples: Some(popgen::TskitSamplesList::Node(&samples[0..x])),
                ..Default::default()
            };
            generate_counts_and_validate(&ts, Some(&options));
        }
    }
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
        let samples = ts.sample_nodes();
        for x in 2..samples.len() - 1 {
            let options = popgen::FromTreeSequenceOptions {
                samples: Some(popgen::TskitSamplesList::Node(&samples[0..x])),
                ..Default::default()
            };
            generate_counts_and_validate(&ts, Some(&options));
        }
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
        let samples = ts.sample_nodes();
        for x in 2..samples.len() - 1 {
            let options = popgen::FromTreeSequenceOptions {
                samples: Some(popgen::TskitSamplesList::Node(&samples[0..x])),
                ..Default::default()
            };
            generate_counts_and_validate(&ts, Some(&options));
        }
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
        let samples = ts.sample_nodes();
        for x in 2..samples.len() - 1 {
            let options = popgen::FromTreeSequenceOptions {
                samples: Some(popgen::TskitSamplesList::Node(&samples[0..x])),
                ..Default::default()
            };
            generate_counts_and_validate(&ts, Some(&options));
        }
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
        let samples = ts.sample_nodes();
        let options = popgen::FromTreeSequenceOptions {
            samples: Some(popgen::TskitSamplesList::Node(samples)),
            ..Default::default()
        };
        generate_counts_and_validate(&ts, Some(&options));

        for x in 2..samples.len() - 1 {
            let options = popgen::FromTreeSequenceOptions {
                samples: Some(popgen::TskitSamplesList::Node(&samples[0..x])),
                ..Default::default()
            };
            generate_counts_and_validate(&ts, Some(&options));
        }
    }
}
