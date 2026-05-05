use crate::geometry::point::Point;

pub struct PointSet {
    points: Vec<Point>,
    dim: usize, // This is added as an enforcement of consistency in a nodeset
}

impl PointSet {
    pub fn new(points: Vec<Point>) -> Result<Self, String> {
        if points.is_empty() {
            return Ok(Self { points, dim: 0 });
        }

        let dim = points[0].coords.len();

        if !points.iter().all(|p| p.coords.len() == dim) {
            return Err("Inconsistent point dimensions".into());
        }

        Ok(Self { points, dim })
    }

    pub fn len(&self) -> usize {
        self.points.len()
    }

    pub fn dim(&self) -> usize {
        self.dim
    }

    pub fn get(&self, i: usize) -> &Point {
        &self.points[i]
    }
}