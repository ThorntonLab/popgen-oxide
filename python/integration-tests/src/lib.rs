use popgen::stats::GlobalStatistic;
use pyo3::prelude::*;

#[pyclass]
struct TreeSequenceHolder {
    ts: tskit::TreeSequence,
}

#[pyclass]
struct SingleSampleCounts {
    counts: popgen::MultiSiteCounts,
}

#[pyclass]
#[repr(transparent)]
struct SingleSampleCountCollection(Vec<SingleSampleCounts>);

#[pymethods]
impl SingleSampleCountCollection {
    pub fn diversity(&self) -> Vec<f64> {
        self.0
            .iter()
            .map(|c| {
                popgen::stats::Diversity::try_from_iter_sites(c.counts.iter())
                    .unwrap()
                    .as_raw()
            })
            .collect::<Vec<_>>()
    }
}

/// A Python module implemented in Rust.
#[pymodule]
mod integration_tests {
    use popgen::stats::GlobalStatistic;
    use pyo3::prelude::*;

    use crate::{SingleSampleCountCollection, SingleSampleCounts, TreeSequenceHolder};

    #[pyfunction]
    fn ts_holder_from_tables(py: Python<'_>, pytables: Py<PyAny>) -> PyResult<TreeSequenceHolder> {
        let shared =
            unsafe { tskit2tskit::SharedTableCollection::new_from_tables(py, pytables).unwrap() };
        let ts = unsafe {
            shared
                .with_tables(|tables| {
                    tables
                        .deepcopy()
                        .unwrap()
                        .tree_sequence(tskit::TreeSequenceFlags::default().build_indexes())
                })
                .unwrap()
        };
        Ok(TreeSequenceHolder { ts })
    }

    #[pyfunction]
    fn counts_from_ts_holder(holder: &TreeSequenceHolder) -> PyResult<SingleSampleCounts> {
        let counts = popgen::MultiSiteCounts::try_from_tree_sequence(
            &holder.ts,
            holder
                .ts
                .node_iter()
                .filter(|n| n.flags().is_sample())
                .map(|n| n.id()),
            None,
        )
        .unwrap();
        Ok(SingleSampleCounts { counts })
    }

    #[pyfunction]
    fn counts_from_ts_holder_single_sample_set(
        holder: &TreeSequenceHolder,
        samples: Vec<i32>,
    ) -> PyResult<SingleSampleCounts> {
        let counts = popgen::MultiSiteCounts::try_from_tree_sequence(
            &holder.ts,
            samples.into_iter().map(|i| i.into()),
            None,
        )
        .unwrap();
        Ok(SingleSampleCounts { counts })
    }

    #[pyfunction]
    /// NOTE: the implementation is extremely inefficient!
    fn counts_from_ts_holder_windowed(
        holder: &TreeSequenceHolder,
        samples: Vec<i32>,
        windows: Vec<f64>,
    ) -> PyResult<SingleSampleCountCollection> {
        let mut vcounts = vec![];
        assert!(!windows.is_empty());
        assert!(windows[0] == 0.0);
        assert!(windows[windows.len() - 1] == holder.ts.tables().sequence_length());
        for w in windows.windows(2) {
            assert!(w[0] < w[1]);
            assert!(w[0] >= 0.);
            assert!(w[1] <= holder.ts.tables().sequence_length());
            let counts = popgen::MultiSiteCounts::try_from_tree_sequence_site_iter(
                &holder.ts,
                samples.iter().map(|i| i.into()),
                holder
                    .ts
                    .site_iter()
                    .filter(|s| s.position() >= w[0] && s.position() < w[1]),
                None,
            )
            .unwrap();
            vcounts.push(SingleSampleCounts { counts });
        }
        Ok(SingleSampleCountCollection(vcounts))
    }

    #[pyfunction]
    fn diversity(counts: &SingleSampleCounts) -> PyResult<f64> {
        Ok(
            popgen::stats::Diversity::try_from_iter_sites(counts.counts.iter())
                .unwrap()
                .as_raw(),
        )
    }
}
