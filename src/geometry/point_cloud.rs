use super::{Set, Space};

// TODO docs
pub struct PointCloud<S: Space> {
    elements: Vec<<S::Set as Set>::Element>,
}

impl<S: Space> PointCloud<S> {
    pub fn new(elements: Vec<<S::Set as Set>::Element>) -> Self {
        Self { elements }
    }

    pub fn len(&self) -> usize {
        self.elements.len()
    }

    pub fn get(&self, i: usize) -> &<S::Set as Set>::Element {
        &self.elements[i]
    }

    pub fn is_empty(&self) -> bool {
        self.elements.is_empty()
    }
}
