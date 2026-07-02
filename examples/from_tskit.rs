use popgen::SampleAlleleCounts;

fn process_nodes_from_ts(ts: &tskit::TreeSequence) {
    let counts = SampleAlleleCounts::try_from_tree_sequence(
        ts,
        ts.node_iter().filter_map(|n| {
            if n.flags().is_sample() {
                Some(n.id())
            } else {
                None
            }
        }),
        None,
    )
    .unwrap();

    counts.iter().for_each(|c| println!("{c:?}"));
}

// Get counts from all nodes found in individuals
fn process_individuals_from_ts(ts: &tskit::TreeSequence) {
    // NOTE: we use an adapter type to take ownersip
    // of individuals as we iterate through them.
    // We implement Iterator on the adapter type
    // to facilitate flattening an iterator over individuals
    // into an iterator over all nodes in individuals.
    //
    // Without the adapter, we would have to copy each
    // slice of nodes into a temporary Vec and convert
    // that into an iterator.
    //
    // Future versions of tskit will contain a feature
    // to do this for you.
    struct IndividualNodeIter<'i> {
        ind: tskit::Individual<'i>,
        current_index: usize,
    }
    impl<'i> Iterator for IndividualNodeIter<'i> {
        type Item = tskit::NodeId;
        fn next(&mut self) -> Option<Self::Item> {
            self.current_index += 1;
            self.ind
                .nodes()
                .and_then(|nodes| nodes.get(self.current_index - 1).copied())
        }
    }
    if ts.individuals().num_rows() > 0 {
        let counts = SampleAlleleCounts::try_from_tree_sequence(
            ts,
            ts.individual_iter().flat_map(|i| IndividualNodeIter {
                ind: i,
                current_index: 0,
            }),
            None,
        )
        .unwrap();

        counts.iter().for_each(|c| println!("{c:?}"));
    } else {
        println!("the individual table is empty...")
    }
}

fn main() {
    let ts = make_tree_sequence();
    process_nodes_from_ts(&ts);
    process_individuals_from_ts(&ts);
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
    process_nodes_from_ts(&ts);
}
