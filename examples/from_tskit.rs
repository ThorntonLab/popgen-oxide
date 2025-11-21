use popgen::MultiSiteCounts;

fn process_ts(ts: &tskit::TreeSequence) {
    let counts = MultiSiteCounts::try_from_tree_sequence(ts, None).unwrap();

    counts.iter().for_each(|c| println!("{c:?}"));
}

fn main() {
    let ts = make_tree_sequence();
    process_ts(&ts);
}

fn make_tree_sequence() -> tskit::TreeSequence {
    // start with tables
    let mut tables = tskit::TableCollection::new(100.0).unwrap();

    // add smaple nodes
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

    // add non-sample nodes
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

    // connect the nodes
    tables.add_edge(0., 100., 4, 0).unwrap();
    tables.add_edge(0., 100., 4, 1).unwrap();
    tables.add_edge(0., 100., 5, 4).unwrap();
    tables.add_edge(0., 100., 5, 2).unwrap();
    tables.add_edge(0., 100., 6, 5).unwrap();
    tables.add_edge(0., 100., 6, 3).unwrap();

    // add some sites and mutations
    let site = tables.add_site(5.0, Some(b"A")).unwrap();
    tables
        .add_mutation(site, 1, tskit::MutationId::NULL, 1.0, Some(b"G"))
        .unwrap();
    let site = tables.add_site(6.0, Some(b"T")).unwrap();
    tables
        .add_mutation(site, 1, tskit::MutationId::NULL, 1.0, Some(b"G"))
        .unwrap();
    tables
        .add_mutation(site, 0, tskit::MutationId::NULL, 1.0, Some(b"C"))
        .unwrap();

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
fn test_process_ts() {
    let ts = make_tree_sequence();
    process_ts(&ts);
}
