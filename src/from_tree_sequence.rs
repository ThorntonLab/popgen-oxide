use crate::{MultiPopulationCounts, MultiSiteCounts, PopgenError, PopgenResult};
use itertools::Itertools;

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
    // NOTE: we manually handle all tskit conventions:
    // * i32 for ids
    // * -1 for "NULL" value
    parent: Vec<i32>,
    num_sample_descendants: Vec<i64>,
    num_mutated_sample_descendants: Vec<i64>,
}

impl TreeData {
    fn new(ts: &tskit::TreeSequence) -> Self {
        Self {
            parent: vec![-1; ts.nodes().num_rows().as_usize()],
            num_sample_descendants: vec![0; ts.nodes().num_rows().as_usize()],
            num_mutated_sample_descendants: vec![0; ts.mutations().num_rows().as_usize()],
        }
    }

    fn update_ancestors(&mut self, p: i32, delta: i64) {
        let mut p = p;
        while p != -1 {
            self.num_sample_descendants[p as usize] += delta;
            debug_assert!(self.num_sample_descendants[p as usize] >= 0);
            p = self.parent[p as usize];
        }
    }

    fn process_output_edge(&mut self, parent: usize, child: usize) {
        debug_assert!(
            self.num_sample_descendants[child] <= self.num_sample_descendants[parent],
            "{parent} ({}) -> {child} ({})",
            self.num_sample_descendants[parent],
            self.num_sample_descendants[child],
        );
        if self.num_sample_descendants[child] > 0 {
            let delta = -self.num_sample_descendants[child];
            self.update_ancestors(self.parent[child], delta);
        }
        self.parent[child] = -1;
    }

    fn process_input_edge(&mut self, parent: usize, child: usize) {
        debug_assert!(self.num_sample_descendants[parent] >= 0);
        debug_assert!(self.num_sample_descendants[child] >= 0);
        self.parent[child] = parent as i32;
        if self.num_sample_descendants[child] > 0 {
            let delta = self.num_sample_descendants[child];
            self.update_ancestors(self.parent[child], delta);
        }
    }

    // This is intended as a truly private fn
    fn get_num_sample_descendants(&mut self, mutation: &tskit::MutationRef<'_>) -> i64 {
        let rv: i64 = self.num_sample_descendants[mutation.node().as_usize()]
            - self.num_mutated_sample_descendants[mutation.id().as_usize()];
        // this is a HARD error representing a serious bug.
        assert!(rv >= 0);
        rv
    }

    #[must_use]
    fn process_mutation<M>(&mut self, mutation: &tskit::MutationRef<'_>, mutation_parent: &M) -> i64
    where
        M: tskit::TableColumn<tskit::MutationId, tskit::MutationId>,
    {
        let nd = self.get_num_sample_descendants(mutation);
        if nd > 0 {
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
        nd
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

fn setup_samples<N>(
    ts: &tskit::TreeSequence,
    samples: N,
) -> Result<(TreeData, i32), crate::PopgenError>
where
    N: Iterator<Item = tskit::NodeId>,
{
    let mut tree_data = TreeData::new(ts);
    let num_nodes = ts.nodes().num_rows().as_usize();
    let num_sampled_genomes = setup_samples_from_node_ids(num_nodes, samples, &mut tree_data)?;
    Ok((tree_data, num_sampled_genomes))
}

fn setup_multi_sample_sets<Outer, Inner>(
    ts: &tskit::TreeSequence,
    samples: Outer,
) -> Result<Vec<(TreeData, i32)>, crate::PopgenError>
where
    Outer: Iterator<Item = Inner>,
    Inner: Iterator<Item = tskit::NodeId>,
{
    let mut rv = vec![];
    for iter in samples {
        let inner = setup_samples(ts, iter)?;
        rv.push(inner)
    }
    Ok(rv)
}

struct SingleSampleSet<'ts> {
    tree_data: TreeData,
    num_sampled_genomes: i64,
    alleles_at_site: Vec<&'ts [u8]>,
    allele_counts: Vec<i64>,
    counts: MultiSiteCounts,
}

struct MultitpleSampleSets<'ts> {
    tree_data: Vec<TreeData>,
    num_sampled_genomes: Vec<i64>,
    alleles_at_site: Vec<&'ts [u8]>,
    allele_counts: Vec<Vec<i64>>,
    num_samples_inheriting_derived_state_at_site: Vec<i64>,
    counts: MultiPopulationCounts,
}

trait SampleSets<'s> {
    type Output: Sized;

    fn process_input_edge(&mut self, parent: usize, child: usize);
    fn process_output_edge(&mut self, parent: usize, child: usize);
    fn initialize_site<'ts, 'a>(
        &'a mut self,
        ts: &'ts tskit::TreeSequence,
        site: tskit::SiteId,
    ) -> PopgenResult<()>
    where
        'ts: 's,
        's: 'a;
    fn process_mutation<'ts, 'a, M>(
        &'a mut self,
        ts: &'ts tskit::TreeSequence,
        mutation_parent: &'a M,
        mutation: tskit::MutationRef<'a>,
    ) -> PopgenResult<()>
    where
        'ts: 's,
        's: 'a,
        M: tskit::TableColumn<tskit::MutationId, tskit::MutationId>;
    fn update_allele_counts(&mut self) -> PopgenResult<()>;
    fn output(self) -> Self::Output;
}

fn setup_alleles_at_site<'ts, 'a>(
    ts: &'ts tskit::TreeSequence,
    site: tskit::SiteId,
    alleles_at_site: &mut Vec<&'a [u8]>,
) -> PopgenResult<()>
where
    'ts: 'a,
{
    alleles_at_site.clear();
    // NOTE: trying to store the derived state
    // from the current_site as a slice runs
    // into lifetime issues because current_site
    // goes away. So what we do instead is get a slice
    // for the same row whose lifetime depends on
    // the tree sequence!
    alleles_at_site.push(*ts.sites().ancestral_state(site).as_ref().ok_or(
        PopgenError::LibraryError("site is missing ancestral state".to_string()),
    )?);
    Ok(())
}

impl<'s> SampleSets<'s> for SingleSampleSet<'s> {
    type Output = MultiSiteCounts;

    fn process_input_edge(&mut self, parent: usize, child: usize) {
        self.tree_data.process_input_edge(parent, child);
    }
    fn process_output_edge(&mut self, parent: usize, child: usize) {
        self.tree_data.process_output_edge(parent, child);
    }

    fn process_mutation<'ts, 'a, M>(
        &'a mut self,
        ts: &'ts tskit::TreeSequence,
        mutation_parent: &'a M,
        mutation: tskit::MutationRef<'a>,
    ) -> PopgenResult<()>
    where
        'ts: 's,
        's: 'a,
        M: tskit::TableColumn<tskit::MutationId, tskit::MutationId>,
    {
        let num_samples_inheriting_derived_state =
            self.tree_data.process_mutation(&mutation, mutation_parent);
        if num_samples_inheriting_derived_state > 0
            && num_samples_inheriting_derived_state < self.num_sampled_genomes
        {
            let derived_state = *ts.mutations().derived_state(mutation.id()).as_ref().ok_or(
                PopgenError::LibraryError("mutation is missing derived state".to_string()),
            )?;
            match self
                .alleles_at_site
                .iter()
                .position(|&x| x == derived_state)
            {
                Some(index) => {
                    if index > 0 {
                        self.allele_counts[index] += num_samples_inheriting_derived_state
                    }
                }
                None => {
                    self.alleles_at_site.push(derived_state);
                    self.allele_counts
                        .push(num_samples_inheriting_derived_state);
                }
            }
        }
        Ok(())
    }

    fn update_allele_counts(&mut self) -> PopgenResult<()> {
        self.allele_counts[0] =
        // TODO: we should simply sum the desired quantity as we go along,
        // eliminating the need for an iteration here.
            (self.num_sampled_genomes) - self.allele_counts.iter().skip(1).sum::<i64>();
        assert!(self.allele_counts[0] >= 0);
        if self
            .allele_counts
            .iter()
            .filter(|&&i| i > 0 && i < self.num_sampled_genomes)
            .count()
            > 1
        {
            self.counts
                .add_site_from_counts(&self.allele_counts, self.num_sampled_genomes as i32)?;
        }
        Ok(())
    }

    fn initialize_site<'ts, 'a>(
        &'a mut self,
        ts: &'ts tskit::TreeSequence,
        site: tskit::SiteId,
    ) -> PopgenResult<()>
    where
        'ts: 's,
        's: 'a,
    {
        setup_alleles_at_site(ts, site, &mut self.alleles_at_site)?;
        self.allele_counts.resize(1, 0);
        Ok(())
    }

    fn output(self) -> Self::Output {
        self.counts
    }
}

impl<'s> SampleSets<'s> for MultitpleSampleSets<'s> {
    type Output = MultiPopulationCounts;

    fn process_input_edge(&mut self, parent: usize, child: usize) {
        self.tree_data
            .iter_mut()
            .for_each(|td| td.process_input_edge(parent, child));
    }
    fn process_output_edge(&mut self, parent: usize, child: usize) {
        self.tree_data
            .iter_mut()
            .for_each(|td| td.process_output_edge(parent, child));
    }

    fn process_mutation<'ts, 'a, M>(
        &'a mut self,
        ts: &'ts tskit::TreeSequence,
        mutation_parent: &'a M,
        mutation: tskit::MutationRef<'a>,
    ) -> PopgenResult<()>
    where
        'ts: 's,
        's: 'a,
        M: tskit::TableColumn<tskit::MutationId, tskit::MutationId>,
    {
        let mut any_sample_sets_polymorphic = false;
        self.num_samples_inheriting_derived_state_at_site.clear();
        for (nd, num_genomes) in self
            .tree_data
            .iter_mut()
            .zip(self.num_sampled_genomes.iter())
            .map(|(tree_data, num_genomes)| {
                (
                    tree_data.process_mutation(&mutation, mutation_parent),
                    num_genomes,
                )
            })
        {
            // Check if mutation is polymorphic in this sample set
            if nd > 0 && nd < *num_genomes {
                any_sample_sets_polymorphic = true;
            }
            self.num_samples_inheriting_derived_state_at_site.push(nd);
        }
        if any_sample_sets_polymorphic {
            let derived_state = *ts.mutations().derived_state(mutation.id()).as_ref().ok_or(
                PopgenError::LibraryError("mutation missing derived state".to_string()),
            )?;
            match self
                .alleles_at_site
                .iter()
                .position(|&x| x == derived_state)
            {
                Some(index) => {
                    if index > 0 {
                        for (i, j) in self
                            .num_samples_inheriting_derived_state_at_site
                            .iter()
                            .enumerate()
                        {
                            self.allele_counts[i][index] += j;
                        }
                    }
                }
                None => {
                    self.alleles_at_site.push(derived_state);
                    for (i, j) in self
                        .num_samples_inheriting_derived_state_at_site
                        .iter()
                        .enumerate()
                    {
                        self.allele_counts[i].push(*j);
                    }
                }
            }
        }
        self.num_samples_inheriting_derived_state_at_site.clear();
        Ok(())
    }

    fn update_allele_counts(&mut self) -> PopgenResult<()> {
        self.num_sampled_genomes
            .iter()
            .enumerate()
            .for_each(|(i, num_sampled_genomes)| {
                self.allele_counts[i][0] =
                    (*num_sampled_genomes) - self.allele_counts[i].iter().skip(1).sum::<i64>()
            });
        assert!(self.allele_counts.iter().all(|v| v[0] >= 0));
        // If ANY of the sample sets are polymorphic,
        // record data for them.
        if self.allele_counts.iter().enumerate().any(|(i, ac)| {
            ac.iter()
                .filter(|&&c| c > 0 && c < self.num_sampled_genomes[i])
                .count()
                > 1
        }) {
            self.counts.extend_populations_from_site(|index| {
                (
                    &self.allele_counts[index],
                    self.num_sampled_genomes[index] as i32,
                )
            })?;
        }
        Ok(())
    }

    fn initialize_site<'ts, 'a>(
        &'a mut self,
        ts: &'ts tskit::TreeSequence,
        site: tskit::SiteId,
    ) -> PopgenResult<()>
    where
        'ts: 's,
        's: 'a,
    {
        setup_alleles_at_site(ts, site, &mut self.alleles_at_site)?;
        self.allele_counts.iter_mut().for_each(|v| v.resize(1, 0));
        Ok(())
    }

    fn output(self) -> Self::Output {
        self.counts
    }
}

fn try_from_tree_sequence_details<'s, S>(
    ts: &'s tskit::TreeSequence,
    parameters: Option<FromTreeSequenceOptions>,
    sample_sets: S,
) -> PopgenResult<S::Output>
where
    S: SampleSets<'s>,
{
    let _parameters = parameters.unwrap_or_default();
    let mut sample_sets = sample_sets;
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
    let mut current_site_index = 0;
    let mut site_iter = ts.site_iter();
    let mut current_site = site_iter.next();
    while i < num_edges && left < ts.tables().sequence_length() {
        while j < num_edges && edges_right[edges_out[j]] == left {
            let edge_parent = edges_parent[edges_out[j]].as_usize();
            let edge_child = edges_child[edges_out[j]].as_usize();
            sample_sets.process_output_edge(edge_parent, edge_child);
            j += 1;
        }
        while i < num_edges && edges_left[edges_in[i]] == left {
            let edge_parent = edges_parent[edges_in[i]].as_usize();
            let edge_child = edges_child[edges_in[i]].as_usize();
            sample_sets.process_input_edge(edge_parent, edge_child);
            i += 1;
        }
        let right = update_right(
            ts.tables().sequence_length().into(),
            i,
            &edges_left,
            &edges_in,
        );
        let right = update_right(right, j, &edges_right, &edges_out);
        while let Some(site_ref) = current_site.as_ref() {
            if site_ref.position() < right {
                sample_sets.initialize_site(ts, site_ref.id())?;

                // NOTE: we process in reverse order because
                // more recent mutations get processed first,
                // allowing the propagation of already-mutated
                // nodes up the tree.
                for mutation in site_ref.mutation_iter().rev() {
                    sample_sets.process_mutation(ts, &mutation_parent, mutation)?;
                }
                sample_sets.update_allele_counts()?;
                current_site = site_iter.next();
            } else {
                break;
            }
        }
        //for current_site in site_iter.peeking_take_while(|site| site.position() < right) {
        //    sample_sets.initialize_site(ts, current_site.id())?;

        //    // NOTE: we process in reverse order because
        //    // more recent mutations get processed first,
        //    // allowing the propagation of already-mutated
        //    // nodes up the tree.
        //    for mutation in current_site.mutation_iter().rev() {
        //        sample_sets.process_mutation(ts, &mutation_parent, mutation)?;
        //    }
        //    sample_sets.update_allele_counts()?;
        //    current_site_index += 1;
        //}
        left = right;
        num_trees += 1;
    }
    // assert_eq!(current_site_index, ts.sites().num_rows().as_usize());
    assert_eq!(num_trees, ts.num_trees());
    Ok(sample_sets.output())
}

pub fn try_from_tree_sequence<N>(
    ts: &tskit::TreeSequence,
    samples: N,
    parameters: Option<FromTreeSequenceOptions>,
) -> PopgenResult<MultiSiteCounts>
where
    N: Iterator<Item = tskit::NodeId>,
{
    let (tree_data, num_sampled_genomes) = setup_samples(ts, samples)?;
    let sample_sets = SingleSampleSet {
        tree_data,
        num_sampled_genomes: num_sampled_genomes as i64,
        alleles_at_site: vec![],
        allele_counts: vec![],
        counts: MultiSiteCounts::default(),
    };
    try_from_tree_sequence_details(ts, parameters, sample_sets)
}

pub fn try_from_tree_sequence_multi<Outer, Inner>(
    ts: &tskit::TreeSequence,
    samples: Outer,
    options: Option<FromTreeSequenceOptions>,
) -> Result<crate::MultiPopulationCounts, PopgenError>
where
    Outer: Iterator<Item = Inner>,
    Inner: Iterator<Item = tskit::NodeId>,
{
    let sample_data = setup_multi_sample_sets(ts, samples)?;
    let counts = MultiPopulationCounts::of_empty_populations(sample_data.len());
    let (tree_data, num_sampled_genomes): (Vec<TreeData>, Vec<i32>) =
        sample_data.into_iter().unzip();
    let num_sampled_genomes: Vec<i64> = num_sampled_genomes.into_iter().map(|i| i as i64).collect();
    let allele_counts = vec![vec![]; tree_data.len()];
    let sample_sets = MultitpleSampleSets {
        tree_data,
        num_sampled_genomes,
        counts,
        allele_counts,
        alleles_at_site: vec![],
        num_samples_inheriting_derived_state_at_site: vec![],
    };
    try_from_tree_sequence_details(ts, options, sample_sets)
}
