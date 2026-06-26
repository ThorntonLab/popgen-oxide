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
    /// If this element exists, the order of elements does not matter.
    /// Because this is a strictly triangular matrix, `x` may not equal `y`.
    /// Panics if this is out of range.
    pub fn get(&self, x: usize, y: usize) -> &T {
        let (x, y) = match x.cmp(&y) {
            Ordering::Greater => (y, x),
            _ => (x, y),
        };
        &self.0[Self::base_for(y) + x]
    }

    /// Extend `self` by providing, for new matrix size `N+1`, the values corresponding to `(i, N+1)` for all `i < N+1`.
    ///
    /// The iterator is always consumed.
    /// The first error produced by the iterator, if any, is propagated.
    /// A rollback to the state before this function was called is not guaranteed.
    ///
    /// If the iterator does not contain precisely `N` elements, the function panics.
    pub fn try_extend<I, E>(&mut self, elems: I) -> Result<(), E>
    where
        I: IntoIterator<Item = Result<T, E>>,
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

#[test]
fn tri_sanity() {
    let mut tri = StrictlyLowerTriangular::<usize>::new();
    tri.try_extend(std::iter::empty::<Result<_, ()>>()).unwrap();
    tri.try_extend(std::iter::once(Ok::<_, ()>(1))).unwrap();

    assert_eq!(tri.get(0, 1), &1);
    assert_eq!(tri.get(1, 0), &1);

    tri.try_extend([2, 3].into_iter().map(Ok::<_, ()>)).unwrap();
    assert_eq!(tri.get(0, 2), &2);
    assert_eq!(tri.get(2, 0), &2);
    assert_eq!(tri.get(1, 2), &3);
    assert_eq!(tri.get(2, 1), &3);
}

#[test]
#[should_panic]
fn tri_panic_1() {
    let mut tri = StrictlyLowerTriangular::<()>::new();
    // should be empty
    tri.try_extend([Ok::<_, ()>(())]).unwrap();
}

#[test]
#[should_panic]
fn tri_panic_2() {
    let mut tri = StrictlyLowerTriangular::<()>::new();
    tri.try_extend(std::iter::empty::<Result<_, ()>>()).unwrap();

    // should have length 1
    tri.try_extend([(), ()].into_iter().map(Ok::<_, ()>))
        .unwrap();
}

#[test]
fn tri_err_propagate() {
    let mut tri = StrictlyLowerTriangular::new();
    tri.try_extend(std::iter::empty::<Result<_, ()>>()).unwrap();
    tri.try_extend(std::iter::once(Ok::<_, ()>(()))).unwrap();

    tri.try_extend([Ok(()), Err(())]).unwrap_err();
}
