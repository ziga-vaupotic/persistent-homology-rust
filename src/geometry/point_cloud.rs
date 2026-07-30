use super::{EuclideanSpace, Space};

// TODO docs
pub struct PointCloud<S: Space> {
    elements: Vec<S::Element>,
}

impl<S: Space> PointCloud<S> {
    pub fn new(elements: Vec<S::Element>) -> Self {
        Self { elements }
    }

    pub fn len(&self) -> usize {
        self.elements.len()
    }

    pub fn get(&self, i: usize) -> &S::Element {
        &self.elements[i]
    }

    pub fn is_empty(&self) -> bool {
        self.elements.is_empty()
    }
}

pub type EuclideanCloud<const N: usize> = PointCloud<EuclideanSpace<N>>;
