use crate::{MultiSiteCounts, PopgenResult};

/// Definition of a single sample set in a [`tskit::TreeSequence`]
pub enum SingleSampleSet<'samples> {
    /// Define samples by **node** identifiers
    Node(&'samples mut dyn Iterator<Item = tskit::NodeId>),
    /// Use all nodes marked as samples
    AllNodes,
    /// Define samples by **individual** identifiers
    Individual(&'samples mut dyn Iterator<Item = tskit::IndividualId>),
    /// Extract nodes from all rows of the individual table
    AllIndividuals,
}

#[derive(Debug, Default)]
/// Options affecting the behavior of
/// [crate::MultiSiteCounts::try_from_tree_sequence]
pub struct FromTreeSequenceOptions {}

fn update_right<P, E>(right: f64, index: usize, position_slice: &P, diff_slice: &E) -> f64
where
    P: tskit::TableColumn<tskit::EdgeId, tskit::Position>,
    E: tskit::TableColumn<tskit::EdgeId, tskit::EdgeId>,
{
    if index < diff_slice.len() {
        let temp = position_slice[diff_slice[index]];
        if temp < right {
            temp.into()
        } else {
            right
        }
    } else {
        right
    }
}

// Poperties of the current tree with respect to
// a given sample set.
// Multi-sample-set tasks will need a Vec of these.
struct TreeData {
    num_sample_descendants: Vec<i64>,
    num_mutated_sample_descendants: Vec<i64>,
}

impl TreeData {
    fn new(ts: &tskit::TreeSequence) -> Self {
        Self {
            num_sample_descendants: vec![0; ts.nodes().num_rows().as_usize()],
            num_mutated_sample_descendants: vec![0; ts.mutations().num_rows().as_usize()],
        }
    }

    fn process_output_edge(&mut self, parent: usize, child: usize) {
        self.num_sample_descendants[parent] -= self.num_sample_descendants[child]
    }

    fn process_input_edge(&mut self, parent: usize, child: usize) {
        self.num_sample_descendants[parent] += self.num_sample_descendants[child]
    }

    // This is intended as a truly private fn
    fn get_num_sample_descendants(&mut self, mutation: &tskit::MutationRef<'_>) -> i64 {
        let rv: i64 = self.num_sample_descendants[mutation.node().as_usize()]
            - self.num_mutated_sample_descendants[mutation.id().as_usize()];
        // this is a HARD error representing a serious bug.
        assert!(rv >= 0);
        rv
    }

    fn process_mutation<M>(
        &mut self,
        mutation: &tskit::MutationRef<'_>,
        mutation_parent: &M,
        alleles_at_site: &[Vec<u8>],
        allele_counts: &mut [i64],
    ) where
        M: tskit::TableColumn<tskit::MutationId, tskit::MutationId>,
    {
        let nd = self.get_num_sample_descendants(mutation);
        if nd > 0 {
            // NOTE: This COULD be an Err, letting the user know they are
            // breaking our expecteations.
            let derived_state = mutation.derived_state().as_ref().unwrap().to_vec();

            match alleles_at_site.iter().position(|x| x == &derived_state) {
                Some(index) if index > 0 => allele_counts[index] += nd,
                Some(_) => (),
                None => unreachable!(),
            };
            // Propagate number of nodes inheriting this mutation up the tree
            let delta = self.num_sample_descendants[mutation.node().as_usize()]
                - self.num_mutated_sample_descendants[mutation.id().as_usize()];
            assert!(!delta.is_negative());
            let mut current_mut_parent = mutation_parent[mutation.id()];
            while !current_mut_parent.is_null() {
                self.num_mutated_sample_descendants[current_mut_parent.as_usize()] += delta;
                current_mut_parent = mutation_parent[current_mut_parent];
            }
        }
    }
}

fn setup_samples_from_node_ids<I>(
    num_nodes: usize,
    iter: I,
    td: &mut TreeData,
) -> Result<i32, crate::PopgenError>
where
    I: Iterator<Item = tskit::NodeId>,
{
    let mut num_sampled_genomes = 0;
    for node_id in iter {
        // Should be an Err condition!
        if node_id == tskit::NodeId::NULL {
            return Err(crate::PopgenError::LibraryError("null node id".to_owned()));
        }
        // Should be an Err condition!
        assert!(node_id.as_usize() < num_nodes);
        if let Some(value) = td.num_sample_descendants.get_mut(node_id.as_usize()) {
            *value += 1;
        } else {
            return Err(crate::PopgenError::LibraryError(format!(
                "node id {node_id} out of range"
            )));
        }
        num_sampled_genomes += 1;
    }
    Ok(num_sampled_genomes)
}

fn setup_samples(
    ts: &tskit::TreeSequence,
    samples: SingleSampleSet<'_>,
) -> Result<(TreeData, i32), crate::PopgenError> {
    let mut tree_data = TreeData::new(ts);
    let num_nodes = ts.nodes().num_rows().as_usize();
    match samples {
        SingleSampleSet::Node(nodes) => {
            let x = setup_samples_from_node_ids(num_nodes, nodes, &mut tree_data)?;
            Ok((tree_data, x))
        }
        SingleSampleSet::AllNodes => {
            let x = setup_samples_from_node_ids(
                num_nodes,
                &mut ts.node_iter().filter_map(|n| {
                    if n.flags().is_sample() {
                        Some(n.id())
                    } else {
                        None
                    }
                }),
                &mut tree_data,
            )?;
            Ok((tree_data, x))
        }
        SingleSampleSet::Individual(individuals) => {
            if ts.individuals().num_rows() == 0 {
                let msg = "tskit::IndividualIds passed for sample list".to_owned()
                    + " but tree sequence has an empty individuals table";
                return Err(crate::PopgenError::LibraryError(msg));
            }
            let nodes = {
                let mut nodes = vec![];
                for individual in individuals {
                    // NOTE: nth is constant-time in tskit > 0.16.3
                    // The next tskit release will also have .individual(id)
                    // for constant-time access.
                    // (not released as of May 5, 2026)
                    if let Some(row) = ts.individual_iter().nth(individual.as_usize()) {
                        if let Some(ind_nodes) = row.nodes() {
                            nodes.extend_from_slice(ind_nodes);
                        } else {
                            // NOTE: is this possible in the tskit data model?
                            return Err(crate::PopgenError::LibraryError(
                                "individual not associated with nodes".to_string(),
                            ));
                        }
                    } else {
                        return Err(crate::PopgenError::LibraryError(format!(
                            "individual id {} out of range",
                            individual
                        )));
                    }
                }
                nodes
            };
            let x = setup_samples_from_node_ids(num_nodes, nodes.into_iter(), &mut tree_data)?;
            Ok((tree_data, x))
        }
        SingleSampleSet::AllIndividuals => {
            if ts.individuals().num_rows() == 0 {
                let msg = "tskit::IndividualIds passed for sample list".to_owned()
                    + " but tree sequence has an empty individuals table";
                return Err(crate::PopgenError::LibraryError(msg));
            }
            let nodes = {
                let mut nodes = vec![];
                for individual in ts.individual_iter() {
                    // NOTE: nth is constant-time in tskit > 0.16.3
                    // The next tskit release will also have .individual(id)
                    // for constant-time access.
                    // (not released as of May 5, 2026)
                    if let Some(ind_nodes) = individual.nodes() {
                        nodes.extend_from_slice(ind_nodes);
                    } else {
                        // NOTE: is this possible in the tskit data model?
                        return Err(crate::PopgenError::LibraryError(
                            "individual not associated with nodes".to_string(),
                        ));
                    }
                }
                nodes
            };
            let x = setup_samples_from_node_ids(num_nodes, nodes.into_iter(), &mut tree_data)?;
            Ok((tree_data, x))
        }
    }
}

pub fn try_from_tree_sequence(
    ts: &tskit::TreeSequence,
    samples: SingleSampleSet<'_>,
    parameters: Option<FromTreeSequenceOptions>,
) -> PopgenResult<MultiSiteCounts> {
    let _parameters = parameters.unwrap_or_default();
    let (mut tree_data, num_sampled_genomes) = setup_samples(ts, samples)?;
    let mut counts = MultiSiteCounts::default();
    let mut left = 0.0;
    let edges_in = ts.edge_insertion_order_column();
    let edges_out = ts.edge_removal_order_column();
    let edges_left = ts.tables().edges().left_column();
    let edges_right = ts.tables().edges().right_column();
    let edges_parent = ts.tables().edges().parent_column();
    let edges_child = ts.tables().edges().child_column();
    let mutation_parent = ts.tables().mutations().parent_column();
    let num_edges = ts.edges().num_rows().as_usize();
    let mut i = 0_usize;
    let mut j = 0_usize;

    let mut num_trees = 0;
    let mut parent = vec![tskit::NodeId::NULL; ts.nodes().num_rows().as_usize()];
    let mut current_site_index = 0;
    let mut alleles_at_site = vec![];
    while i < num_edges && left < ts.tables().sequence_length() {
        while j < num_edges && edges_right[edges_out[j]] == left {
            let edge_parent = edges_parent[edges_out[j]].as_usize();
            let edge_child = edges_child[edges_out[j]].as_usize();
            tree_data.process_output_edge(edge_parent, edge_child);
            parent[edges_child[edges_out[j]].as_usize()] = tskit::NodeId::NULL;
            j += 1;
        }
        while i < num_edges && edges_left[edges_in[i]] == left {
            parent[edges_child[edges_in[i]].as_usize()] = edges_parent[edges_in[i]];
            let edge_parent = edges_parent[edges_in[i]].as_usize();
            let edge_child = edges_child[edges_in[i]].as_usize();
            tree_data.process_input_edge(edge_parent, edge_child);
            i += 1;
        }
        let right = update_right(
            ts.tables().sequence_length().into(),
            i,
            &edges_left,
            &edges_in,
        );
        let right = update_right(right, j, &edges_right, &edges_out);
        for current_site in ts
            .site_iter()
            .skip(current_site_index)
            .take_while(|site| site.position() < right)
        {
            alleles_at_site.clear();
            // Hard error intentional -- these calcs cannot be done w/o state data
            // TODO: this missing state might mean something (e.g. insertion); figure this out later
            alleles_at_site.push(current_site.ancestral_state().as_ref().unwrap().to_vec());

            // Extract derived states
            alleles_at_site.extend(
                current_site
                    .mutation_iter()
                    .rev()
                    .map(|m| m.derived_state().unwrap().to_vec()),
            );

            let mut allele_counts = vec![0_i64; alleles_at_site.len()];

            // NOTE: we process in reverse order because
            // more recent mutations get processed first,
            // allowing the propagation of already-mutated
            // nodes up the tree.
            for mutation in current_site.mutation_iter().rev() {
                tree_data.process_mutation(
                    &mutation,
                    &mutation_parent,
                    &alleles_at_site,
                    &mut allele_counts,
                );
            }
            allele_counts[0] =
                (num_sampled_genomes as i64) - allele_counts.iter().skip(1).sum::<i64>();
            assert!(allele_counts[0] >= 0);
            if allele_counts
                .iter()
                .filter(|&&i| i > 0 && i < num_sampled_genomes as i64)
                .count()
                > 1
            {
                // Filter out any alleles monomorphic in our sample
                // NOTE: in general, we'd only do this for alleles
                // monomorphic across ALL sample sets
                // TODO: why are we keeping the ancestral allele
                // even if it has counts of 0?
                // * skipping the ancestral allele when it has
                //   a count of 0 causes tests to fail.
                // * From the VCF data, the first element is the
                //   referenence allele? Is that the issue?
                let allele_counts = {
                    let mut temp = vec![allele_counts[0]];
                    temp.extend(
                        allele_counts
                            .iter()
                            .skip(1)
                            .cloned()
                            .filter(|&i| i > 0 && i < num_sampled_genomes as i64),
                    );
                    temp
                };
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
    assert_eq!(current_site_index, ts.sites().num_rows().as_usize());
    assert_eq!(num_trees, ts.num_trees());
    Ok(counts)
}
