use std::cmp::Ordering;

#[derive(Ord, PartialOrd, Eq, PartialEq, Copy, Clone, Debug, Hash, Default)]
pub struct UnorderedPair<T: Ord>(pub T, pub T);

impl<T: Ord> UnorderedPair<T> {
    pub fn new(a: T, b: T) -> Self {
        match &a.cmp(&b) {
            Ordering::Less => Self(a, b),
            // stable
            Ordering::Equal => Self(a, b),
            Ordering::Greater => Self(b, a),
        }
    }
}

impl<T: Ord> From<(T, T)> for UnorderedPair<T> {
    fn from((a, b): (T, T)) -> Self {
        Self::new(a, b)
    }
}
