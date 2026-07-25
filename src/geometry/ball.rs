use crate::geometry::Point;

#[derive(Clone)]
pub struct Ball {
    // closed ie. with boundary
    pub center: Point,
    pub radius: f64,
}

impl Ball {
    pub fn new(center: Point, radius: f64) -> Self {
        Self { center, radius }
    }

    pub fn o(&self) -> &Point {
        &self.center
    }

    pub fn r(&self) -> f64 {
        self.radius
    }
}
