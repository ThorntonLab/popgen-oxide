use crate::{MultiSiteCounts, PopgenResult};
use tskit::MutationId;

// NOTES:
// 1. tskit could use Add for all id types!
//    Oh, they do exist but they are for Self.
//    Should refactor to make them accept Into<Self>
// 2. rust's denial of integer casts can make this kind of "data-driven"
//    API annoying if they come from other languages.
//    tskit defines "size_type" as u32 on 64 bit systems
//    and on 32 bit systems unless you compile with "big tables"
//    enabled. (tskit-rust does not support "big tables" right now)
//    This is all "fine and dandy" because C allows implicit casts for any and
//    all integer types so that x[i] makes you not go insane
//    (unless you promote the cast warning to an error).
//    In rust, dealing w/3rd party libs like this has
//    poor ergonomics

#[derive(Default)]
struct FromTreeSequenceOptions {}

fn update_right(
    right: f64,
    index: u64,
    position_slice: &[tskit::Position],
    diff_slice: &[tskit::EdgeId],
) -> f64 {
    let index = index as usize;
    if index < diff_slice.len() {
        let temp = position_slice[diff_slice[index].as_usize()];
        if temp < right {
            temp.into()
        } else {
            right
        }
    } else {
        right
    }
}

fn try_from_tree_sequence(
    ts: &tskit::TreeSequence,
    _parameters: Option<FromTreeSequenceOptions>,
) -> PopgenResult<MultiSiteCounts> {
    let mut counts = MultiSiteCounts::default();
    let mut left = 0.0;
    let mut right = f64::from(ts.tables().sequence_length());
    // NOTE: we need TreeSequence to be able to provide these
    // indexes w/o going thru the Option b/c you CANNOT make a
    // ts from unindexed tables!!
    // the absence of edge orderings unwrapped below (because the tables aren't indexed)
    // is a data model error
    let edges_in = ts.tables().edge_insertion_order().unwrap();
    let edges_out = ts.tables().edge_removal_order().unwrap();
    let edges_left = ts.tables().edges().left_slice();
    let edges_right = ts.tables().edges().right_slice();
    let edges_parent = ts.tables().edges().parent_slice();
    let edges_child = ts.tables().edges().child_slice();
    let site_pos = ts.tables().sites().position_slice();
    let mutation_site = ts.tables().mutations().site_slice();
    let mutation_node = ts.tables().mutations().node_slice();
    let mutation_parent = ts.tables().mutations().parent_slice();
    let num_edges = ts.edges().num_rows();
    let mut i = 0;
    let mut j = 0;

    let mut num_trees = 0;
    let mut num_sample_descendants = vec![0_i64; ts.nodes().num_rows().as_usize()];
    let mut num_mutated_sample_descendants = vec![0_i64; ts.mutations().num_rows().as_usize()];
    let mut parent = vec![tskit::NodeId::NULL; ts.nodes().num_rows().as_usize()];
    for s in ts.nodes().iter().filter(|i| i.flags.is_sample()) {
        num_sample_descendants[s.id.as_usize()] = 1;
    }
    let mut current_site_index = 0;
    let mut current_mutation_index = 0;
    let mut alleles_at_site = vec![];
    while i < num_edges && left < ts.tables().sequence_length() {
        while j < num_edges && edges_right[edges_out[j as usize].as_usize()] == left {
            num_sample_descendants[edges_parent[edges_out[j as usize].as_usize()].as_usize()] -=
                num_sample_descendants[edges_child[edges_out[j as usize].as_usize()].as_usize()];
            parent[edges_child[edges_out[j as usize].as_usize()].as_usize()] = tskit::NodeId::NULL;
            j += 1;
        }
        while i < num_edges && edges_left[edges_in[i as usize].as_usize()] == left {
            parent[edges_child[edges_in[i as usize].as_usize()].as_usize()] =
                edges_parent[edges_in[i as usize].as_usize()];
            num_sample_descendants[edges_parent[edges_in[i as usize].as_usize()].as_usize()] +=
                num_sample_descendants[edges_child[edges_in[i as usize].as_usize()].as_usize()];
            i += 1;
        }
        while current_site_index < ts.sites().num_rows()
            && site_pos[current_site_index as usize] < right
        {
            alleles_at_site.clear();
            alleles_at_site.push(
                ts.sites()
                    .ancestral_state(current_site_index as i32)
                    // Hard error intentional -- these calcs cannot be done w/o state data
                    // TODO: this missing state might mean something (e.g. insertion); figure this out later
                    .unwrap(),
            );
            let mut allele_counts = vec![0_i64];
            while current_mutation_index < ts.mutations().num_rows()
                // Dang, tskit integer types can get frustrating
                && mutation_site[current_mutation_index as usize] == current_site_index as i32
            {
                let temp = mutation_site[current_mutation_index as usize..]
                    .iter()
                    .take_while(|&&site| site == (current_site_index as i32))
                    .count();
                // Experimental code follows
                let mut mnode = None;
                for mutation_index in
                    (current_mutation_index..current_mutation_index + (temp as u64)).rev()
                {
                    if let Some(mut_node) = mnode {
                        if mutation_node[mutation_index as usize] != mut_node {
                            mnode = None;
                        }
                    }

                    if mnode.is_none() {
                        let current_mut_node = mutation_node[mutation_index as usize];
                        let nd = num_sample_descendants[current_mut_node.as_usize()]
                            .checked_sub(num_mutated_sample_descendants[mutation_index as usize])
                            // again -- this is a HARD error representing a serious bug.
                            .unwrap();
                        assert!(nd >= 0, "nd = {nd} at {current_mut_node:?}");
                        if nd > 0 {
                            let derived_state = ts
                                .mutations()
                                .derived_state(mutation_index as i32)
                                // Hard error intentional -- these calcs cannot be done w/o state data
                                // TODO: this might mean something, but out of scope for now
                                .unwrap();

                            if let Some(index) =
                                alleles_at_site.iter().position(|&x| x == derived_state)
                            {
                                if index > 0 {
                                    // NOT the ancestral state!
                                    allele_counts[index] += nd
                                }
                            } else {
                                alleles_at_site.push(derived_state);
                                allele_counts.push(nd)
                            }
                            let delta = num_sample_descendants[current_mut_node.as_usize()]
                                - num_mutated_sample_descendants[mutation_index as usize];
                            assert!(!delta.is_negative());
                            let mut current_mut_parent = mutation_parent[mutation_index as usize];
                            while !current_mut_parent.is_null() {
                                num_mutated_sample_descendants[current_mut_parent.as_usize()] +=
                                    delta;
                                current_mut_parent = mutation_parent[current_mut_parent.as_usize()];
                            }
                        }
                        mnode = Some(current_mut_node);
                    }
                }
                current_mutation_index += temp as u64;
            }
            // Easier way?
            allele_counts[0] =
                (u64::from(ts.num_samples()) as i64) - allele_counts.iter().skip(1).sum::<i64>();
            assert!(allele_counts[0] >= 0);
            if allele_counts
                .iter()
                .filter(|&&i| i > 0 && (i as u64) < ts.num_samples())
                .count()
                > 1
            {
                counts.add_site_from_counts(&allele_counts, allele_counts.len() as i32);
            }
            current_site_index += 1;
        }
        right = update_right(
            ts.tables().sequence_length().into(),
            i,
            edges_left,
            edges_in,
        );
        right = update_right(right, j, edges_right, edges_out);
        left = right;
        num_trees += 1;
    }
    assert_eq!(current_site_index, ts.sites().num_rows());
    assert_eq!(current_mutation_index, ts.mutations().num_rows());
    assert_eq!(num_trees, ts.num_trees());
    Ok(counts)
}

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
                    MutationId::NULL,
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
    let counts = try_from_tree_sequence(&ts, None).unwrap();
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
    let counts = try_from_tree_sequence(&ts, None).unwrap();
    assert_eq!(counts.iter().count(), 1);
    counts
        .iter()
        .take(1)
        .for_each(|c| assert_eq!(c.counts, &[1, 1]));
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
    let counts = try_from_tree_sequence(&ts, None).unwrap();
    assert_eq!(counts.iter().count(), 2);
    counts
        .iter()
        .take(1)
        .for_each(|c| assert_eq!(c.counts, &[1, 1]));
    counts
        .iter()
        .skip(1)
        .take(1)
        .for_each(|c| assert_eq!(c.counts, &[0, 1, 1]));
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
    let counts = try_from_tree_sequence(&ts, None).unwrap();
    assert_eq!(counts.iter().count(), 1);
    counts.iter().for_each(|c| assert_eq!(c.counts, &[2, 2]));
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
    let counts = try_from_tree_sequence(&ts, None).unwrap();
    assert_eq!(counts.iter().count(), 1);
    counts.iter().for_each(|c| assert_eq!(c.counts, &[2, 1, 1]));
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
    let counts = try_from_tree_sequence(&ts, None).unwrap();
    assert_eq!(counts.iter().count(), 2);
    assert_eq!(counts.get(0).unwrap().counts, &[2, 1, 1]);
    assert_eq!(counts.get(1).unwrap().counts, &[2, 1, 1]);
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
    let counts = try_from_tree_sequence(&ts, None).unwrap();
    assert_eq!(counts.iter().count(), 2);
    assert_eq!(counts.get(0).unwrap().counts, &[1, 3]);
    assert_eq!(counts.get(1).unwrap().counts, &[2, 1, 1]);
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
    let counts = try_from_tree_sequence(&ts, None).unwrap();
    assert_eq!(counts.iter().count(), 2);
    assert_eq!(counts.get(0).unwrap().counts, &[1, 3]);
    assert_eq!(counts.get(1).unwrap().counts, &[3, 1]);
}

#[test]
fn test_8() {
    let site0 = SiteData::new(60.0, "G", vec![]);
    let ts = make_test_data(make_two_different_four_sample_trees, vec![site0]);
    let counts = try_from_tree_sequence(&ts, None).unwrap();
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
    let counts = try_from_tree_sequence(&ts, None).unwrap();
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
    let _ = try_from_tree_sequence(&ts, None).unwrap();
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
    let _ = try_from_tree_sequence(&ts, None).unwrap();
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
    let counts = try_from_tree_sequence(&ts, None).unwrap();
    assert_eq!(counts.iter().count(), 1);
    counts.iter().for_each(|c| assert_eq!(c.counts, &[1, 2, 1]));
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
    let counts = try_from_tree_sequence(&ts, None).unwrap();
    assert_eq!(counts.iter().count(), 1);
    counts
        .iter()
        .for_each(|c| assert_eq!(c.counts, &[1, 1, 1, 1, 1, 1]));
}
