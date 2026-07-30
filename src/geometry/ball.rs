use super::{EuclideanSpace, MetricSpace};

// TODO docs
#[derive(Clone)]
pub struct Ball<M: MetricSpace> {
    pub center: M::Element,
    pub radius: f64,
}

impl<M: MetricSpace> Ball<M> {
    pub fn new(center: M::Element, radius: f64) -> Self {
        Self { center, radius }
    }

    pub fn o(&self) -> &M::Element {
        &self.center
    }

    pub fn r(&self) -> f64 {
        self.radius
    }

    pub fn contains(&self, a: &M::Element) -> bool {
        M::distance(self.o(), a) <= self.radius
    }
}

pub type EuclideanBall<const N: usize> = Ball<EuclideanSpace<N>>;
