use varistat::stats::StatRepresentation;
use varistat::stats::UnpolarisedSiteStat;
use pyo3::prelude::*;

#[pyclass]
struct TreeSequenceHolder {
    ts: tskit::TreeSequence,
}

#[pyclass]
struct SingleSampleCounts {
    counts: varistat::SampleAlleleCounts,
}

#[pyclass]
#[repr(transparent)]
struct SingleSampleCountCollection(Vec<SingleSampleCounts>);

#[pymethods]
impl SingleSampleCountCollection {
    /// For the purposes of testing, we treat empty count objects
    /// as having a diversity of 0.0
    pub fn diversity(&self) -> Vec<f64> {
        self.0
            .iter()
            .map(
                |c| match varistat::stats::Diversity::try_from_iter_sites(c.counts.iter()) {
                    Ok(div) => div.as_raw(),
                    Err(varistat::VaristatError::EmptySiteCounts) => {
                        varistat::stats::Diversity::default().as_raw()
                    }
                    Err(e) => panic!("unexpected error {e:?}"),
                },
            )
            .collect::<Vec<_>>()
    }
}

/// A Python module implemented in Rust.
#[pymodule]
mod integration_tests {
    use varistat::stats::StatRepresentation;
    use varistat::stats::UnpolarisedSiteStat;
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
        let counts = varistat::SampleAlleleCounts::try_from_tree_sequence(
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
        let counts = varistat::SampleAlleleCounts::try_from_tree_sequence(
            &holder.ts,
            samples.into_iter().map(|i| i.into()),
            None,
        )
        .unwrap();
        Ok(SingleSampleCounts { counts })
    }

    /// This is windows the "tskit-python" way,
    /// meaning that the windows must span the entire sequence length
    /// of the input
    #[pyfunction]
    fn counts_from_ts_holder_windowed(
        holder: &TreeSequenceHolder,
        samples: Vec<i32>,
        windows: Vec<f64>,
    ) -> PyResult<SingleSampleCountCollection> {
        assert!(!windows.is_empty());
        assert!(windows[0] == 0.0);
        assert!(windows[windows.len() - 1] == holder.ts.tables().sequence_length());
        let counts = varistat::SampleAlleleCounts::try_from_tree_sequence_windows(
            &holder.ts,
            samples.iter().map(|i| i.into()),
            windows.windows(2).map(|w| (w[0], w[1])),
            None,
        )
        .unwrap();
        let vcounts = counts
            .into_iter()
            .map(|c| SingleSampleCounts { counts: c })
            .collect::<Vec<_>>();
        Ok(SingleSampleCountCollection(vcounts))
    }

    /// More general windowing fn.
    #[pyfunction]
    fn counts_from_ts_holder_windowed_general(
        holder: &TreeSequenceHolder,
        samples: Vec<i32>,
        windows: Vec<(f64, f64)>,
    ) -> PyResult<SingleSampleCountCollection> {
        assert!(!windows.is_empty());
        let counts = varistat::SampleAlleleCounts::try_from_tree_sequence_windows(
            &holder.ts,
            samples.iter().map(|i| i.into()),
            windows.into_iter(),
            None,
        )
        .unwrap();
        let vcounts = counts
            .into_iter()
            .map(|c| SingleSampleCounts { counts: c })
            .collect::<Vec<_>>();
        Ok(SingleSampleCountCollection(vcounts))
    }

    #[pyfunction]
    /// For the purposes of testing, we treat empty count objects
    /// as having a diversity of 0.0
    fn diversity(counts: &SingleSampleCounts) -> PyResult<f64> {
        let div = match varistat::stats::Diversity::try_from_iter_sites(counts.counts.iter()) {
            Ok(div) => div.as_raw(),
            Err(varistat::PopgenError::EmptySiteCounts) => {
                varistat::stats::Diversity::default().as_raw()
            }
            Err(e) => panic!("unexpected error {e:?}"),
        };
        Ok(div)
    }
}
