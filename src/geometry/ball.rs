use super::{MetricSpace, Set};

// TODO docs
#[derive(Clone)]
pub struct Ball<M: MetricSpace> {
    pub center: <M::Set as Set>::Element,
    pub radius: f64,
}

impl<M: MetricSpace> Ball<M> {
    pub fn new(center: <M::Set as Set>::Element, radius: f64) -> Self {
        Self { center, radius }
    }

    pub fn o(&self) -> &<M::Set as Set>::Element {
        &self.center
    }

    pub fn r(&self) -> f64 {
        self.radius
    }

    pub fn contains(&self, a: &<M::Set as Set>::Element) -> bool {
        M::distance(self.o(), a) <= self.radius
    }
}
