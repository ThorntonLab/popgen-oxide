use std::cmp::Ordering;

/// A matrix of size `N`x`N` for which any `(i, j)` is empty if `i <= j`.
#[derive(Debug, Clone)]
pub(crate) struct StrictlyLowerTriangular<T>(Vec<T>, usize);

impl<T> StrictlyLowerTriangular<T> {
    pub fn new() -> Self {
        Self(Vec::new(), 0)
    }

    fn base_for(row: usize) -> usize {
        // assume wrapping is the default; because we multiply by n, this is still correct for n=0
        (row * (row.wrapping_sub(1))) / 2
    }

    /// Get the element corresponding to the unordered pair (x, y).
    /// Panics if this is out of range.
    pub fn get(&self, x: usize, y: usize) -> &T {
        let (x, y) = match x.cmp(&y) {
            Ordering::Greater => (y, x),
            _ => (x, y)
        };
        &self.0[Self::base_for(y) + x]
    }

    /// Extend `self` by providing, for new matrix size `N+1`, the values corresponding to `(i, N+1)` for all `i < N+1`.
    /// If the iterator does not contain precisely `N` elements, the function panics.
    /// The iterator is always consumed. The first error produced by the iterator, if any, is propagated.
    pub fn try_extend<I, E>(&mut self, elems: I) -> Result<(), E>
    where
        I: Iterator<Item = Result<T, E>>,
    {
        self.0.reserve(self.1);
        let expected_new_len = self.0.len() + self.1;
        for e in elems {
            self.0.push(e?);
        }

        if self.0.len() != expected_new_len {
            panic!("StrictlyLowerTriangular::try_extend: iterator wrong length");
        }

        self.1 += 1;
        Ok(())
    }
}
