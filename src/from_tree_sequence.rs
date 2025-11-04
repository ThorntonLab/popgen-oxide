use crate::{MultiSiteCounts, PopgenResult};

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
/// Options affecting the behavior of
/// [crate::MultiSiteCounts::try_from_tree_sequence]
pub struct FromTreeSequenceOptions {}

fn update_right(
    right: f64,
    index: usize,
    position_slice: &[tskit::Position],
    diff_slice: &[tskit::EdgeId],
) -> f64 {
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

pub fn try_from_tree_sequence(
    ts: &tskit::TreeSequence,
    _parameters: Option<FromTreeSequenceOptions>,
) -> PopgenResult<MultiSiteCounts> {
    let mut counts = MultiSiteCounts::default();
    let mut left = 0.0;
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
    let num_edges = ts.edges().num_rows().as_usize();
    let mut i = 0_usize;
    let mut j = 0_usize;

    let mut num_trees = 0;
    let mut num_sample_descendants = vec![0_i64; ts.nodes().num_rows().as_usize()];
    let mut num_mutated_sample_descendants = vec![0_i64; ts.mutations().num_rows().as_usize()];
    let mut parent = vec![tskit::NodeId::NULL; ts.nodes().num_rows().as_usize()];
    let mut num_sampled_genomes = 0_i32;
    for s in ts.nodes().iter().filter(|i| i.flags.is_sample()) {
        num_sample_descendants[s.id.as_usize()] = 1;
        num_sampled_genomes += 1;
    }
    let mut current_site_index = 0;
    let mut current_mutation_index = 0;
    let mut alleles_at_site = vec![];
    while i < num_edges && left < ts.tables().sequence_length() {
        while j < num_edges && edges_right[edges_out[j].as_usize()] == left {
            num_sample_descendants[edges_parent[edges_out[j].as_usize()].as_usize()] -=
                num_sample_descendants[edges_child[edges_out[j].as_usize()].as_usize()];
            parent[edges_child[edges_out[j].as_usize()].as_usize()] = tskit::NodeId::NULL;
            j += 1;
        }
        while i < num_edges && edges_left[edges_in[i].as_usize()] == left {
            parent[edges_child[edges_in[i].as_usize()].as_usize()] =
                edges_parent[edges_in[i].as_usize()];
            num_sample_descendants[edges_parent[edges_in[i].as_usize()].as_usize()] +=
                num_sample_descendants[edges_child[edges_in[i].as_usize()].as_usize()];
            i += 1;
        }
        let right = update_right(
            ts.tables().sequence_length().into(),
            i,
            edges_left,
            edges_in,
        );
        let right = update_right(right, j, edges_right, edges_out);
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
                // this won't panic because our counts are ultimately derived from a collection of
                // alleles, which always obeys the required properties
                counts
                    .add_site_from_counts(&allele_counts, num_sampled_genomes)
                    .unwrap();
            }
            current_site_index += 1;
        }
        left = right;
        num_trees += 1;
    }
    assert_eq!(current_site_index, ts.sites().num_rows());
    assert_eq!(current_mutation_index, ts.mutations().num_rows());
    assert_eq!(num_trees, ts.num_trees());
    Ok(counts)
}
