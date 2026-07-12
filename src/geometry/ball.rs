

use crate::geometry::Point;



pub struct Ball { // closed ie. with boundary
    pub center : Point,
    pub radius : f64
}


impl Ball {

    pub fn new(center : Point, radius : f64) -> Self {
        Self { center : center, radius : radius }
    }


    pub fn o(&self) -> &Point {
        &self.center
    }


    pub fn r(&self) -> f64 {
        self.radius
    }


    pub fn contains(&self, point : &Point) -> bool {
        self.center.dim() > 0 && self.center.distance(point) <= self.radius
    }

}
