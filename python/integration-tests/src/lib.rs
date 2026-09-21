use popgen::stats::StatRepresentation;
use popgen::stats::UnpolarisedSiteStat;
use pyo3::prelude::*;

#[pyclass]
struct TreeSequenceHolder {
    ts: tskit::TreeSequence,
}

#[pyclass]
struct SingleSampleCounts {
    counts: popgen::SampleAlleleCounts,
}

#[pyclass]
struct Fstatistics {
    fstats: popgen::stats::FStatistics,
}

#[pymethods]
impl Fstatistics {
    pub fn f2(&self, set1: usize, set2: usize) -> PyResult<f64> {
        let f2 = self.fstats.f2(set1, set2).unwrap();
        Ok(f2)
    }

    pub fn diversity(&self, set: usize) -> PyResult<f64> {
        Ok(self.fstats.pi_within(set).unwrap())
    }
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
            .map(|c| {
                match popgen::stats::Diversity::try_from_iter_sites(
                    c.counts.iter_sample_set(0).unwrap(),
                ) {
                    Ok(div) => div.as_raw(),
                    Err(popgen::PopgenError::EmptySiteCounts) => {
                        popgen::stats::Diversity::default().as_raw()
                    }
                    Err(e) => panic!("unexpected error {e:?}"),
                }
            })
            .collect::<Vec<_>>()
    }
}

/// A Python module implemented in Rust.
#[pymodule]
mod integration_tests {
    use popgen::stats::StatRepresentation;
    use popgen::stats::UnpolarisedSiteStat;
    use pyo3::prelude::*;

    use crate::{Fstatistics, SingleSampleCountCollection, SingleSampleCounts, TreeSequenceHolder};

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
        let counts = popgen::SampleAlleleCounts::try_from_tree_sequence(
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
        let counts = popgen::SampleAlleleCounts::try_from_tree_sequence(
            &holder.ts,
            samples.into_iter().map(|i| i.into()),
            None,
        )
        .unwrap();
        Ok(SingleSampleCounts { counts })
    }

    #[pyfunction]
    fn counts_from_ts_holder_multi_sample_sets(
        holder: &TreeSequenceHolder,
        sample_sets: Vec<Vec<i32>>,
    ) -> PyResult<SingleSampleCounts> {
        let counts = popgen::SampleAlleleCounts::try_multi_sample_set_from_tree_sequence(
            &holder.ts,
            sample_sets
                .into_iter()
                .map(|i| i.into_iter().map(|j| j.into())),
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
        let counts = popgen::SampleAlleleCounts::try_from_tree_sequence_windows(
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
        let counts = popgen::SampleAlleleCounts::try_from_tree_sequence_windows(
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
        let div = match popgen::stats::Diversity::try_from_iter_sites(
            counts.counts.iter_sample_set(0).unwrap(),
        ) {
            Ok(div) => div.as_raw(),
            Err(popgen::PopgenError::EmptySiteCounts) => {
                popgen::stats::Diversity::default().as_raw()
            }
            Err(e) => panic!("unexpected error {e:?}"),
        };
        Ok(div)
    }

    #[pyfunction]
    fn fstats(counts: &SingleSampleCounts) -> PyResult<Fstatistics> {
        let fstats =
            popgen::stats::FStatistics::try_from_sample_sets(&counts.counts, |_| Some(1.)).unwrap();
        Ok(Fstatistics { fstats })
    }
}
