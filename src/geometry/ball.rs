use crate::geometry::Point;

/// A closed ball in point space.
///
/// The ball is defined by a center point and a radius.
#[derive(Clone)]
pub struct Ball {
    // closed ie. with boundary
    pub center: Point,
    pub radius: f64,
}

impl Ball {
    /// Create a new closed ball with the given center and radius.
    pub fn new(center: Point, radius: f64) -> Self {
        Self { center, radius }
    }

    /// Return the center of the ball.
    pub fn o(&self) -> &Point {
        &self.center
    }

    /// Return the radius of the ball.
    pub fn r(&self) -> f64 {
        self.radius
    }
}
