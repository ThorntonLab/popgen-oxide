// NOTE: these tests require compiling
// with the tskit feature

use popgen::tskit;

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

#[test]
fn test_0() {
    let ts = make_test_data(make_two_sample_tree, vec![]);
    let counts = popgen::MultiSiteCounts::try_from_tree_sequence(&ts, None).unwrap();
    assert_eq!(counts.iter().count(), 0)
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
    let counts = popgen::MultiSiteCounts::try_from_tree_sequence(&ts, None).unwrap();
    assert_eq!(counts.iter().count(), 1);
    counts
        .iter()
        .take(1)
        .for_each(|c| assert_eq!(c.counts(), &[1, 1]));
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
    assert_eq!(ts.tables().sites().num_rows(), 2);
    assert_eq!(ts.tables().mutations().num_rows(), 3);
    let counts = popgen::MultiSiteCounts::try_from_tree_sequence(&ts, None).unwrap();
    assert_eq!(counts.iter().count(), 2);
    counts
        .iter()
        .take(1)
        .for_each(|c| assert_eq!(c.counts(), &[1, 1]));
    counts
        .iter()
        .skip(1)
        .take(1)
        .for_each(|c| assert_eq!(c.counts(), &[0, 1, 1]));
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
    let counts = popgen::MultiSiteCounts::try_from_tree_sequence(&ts, None).unwrap();
    assert_eq!(counts.iter().count(), 1);
    counts.iter().for_each(|c| assert_eq!(c.counts(), &[2, 2]));
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
    let counts = popgen::MultiSiteCounts::try_from_tree_sequence(&ts, None).unwrap();
    assert_eq!(counts.iter().count(), 1);
    counts
        .iter()
        .for_each(|c| assert_eq!(c.counts(), &[2, 1, 1]));
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
    let counts = popgen::MultiSiteCounts::try_from_tree_sequence(&ts, None).unwrap();
    assert_eq!(counts.iter().count(), 2);
    assert_eq!(counts.get(0).unwrap().counts(), &[2, 1, 1]);
    assert_eq!(counts.get(1).unwrap().counts(), &[2, 1, 1]);
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
    let counts = popgen::MultiSiteCounts::try_from_tree_sequence(&ts, None).unwrap();
    assert_eq!(counts.iter().count(), 2);
    assert_eq!(counts.get(0).unwrap().counts(), &[1, 3]);
    assert_eq!(counts.get(1).unwrap().counts(), &[2, 1, 1]);
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
    let counts = popgen::MultiSiteCounts::try_from_tree_sequence(&ts, None).unwrap();
    assert_eq!(counts.iter().count(), 2);
    assert_eq!(counts.get(0).unwrap().counts(), &[1, 3]);
    assert_eq!(counts.get(1).unwrap().counts(), &[3, 1]);
}

#[test]
fn test_8() {
    let site0 = SiteData::new(60.0, "G", vec![]);
    let ts = make_test_data(make_two_different_four_sample_trees, vec![site0]);
    let counts = popgen::MultiSiteCounts::try_from_tree_sequence(&ts, None).unwrap();
    assert_eq!(counts.iter().count(), 0);
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
    let counts = popgen::MultiSiteCounts::try_from_tree_sequence(&ts, None).unwrap();
    assert_eq!(counts.iter().count(), 0);
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
    let _ = popgen::MultiSiteCounts::try_from_tree_sequence(&ts, None).unwrap();
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
    let _ = popgen::MultiSiteCounts::try_from_tree_sequence(&ts, None).unwrap();
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
    let counts = popgen::MultiSiteCounts::try_from_tree_sequence(&ts, None).unwrap();
    assert_eq!(counts.iter().count(), 1);
    counts
        .iter()
        .for_each(|c| assert_eq!(c.counts(), &[1, 2, 1]));
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
    let counts = popgen::MultiSiteCounts::try_from_tree_sequence(&ts, None).unwrap();
    assert_eq!(counts.iter().count(), 1);
    counts
        .iter()
        .for_each(|c| assert_eq!(c.counts(), &[1, 1, 1, 1, 1, 1]));
}
