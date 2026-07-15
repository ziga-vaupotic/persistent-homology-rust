


use crate::geometry::Point;

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


    pub fn new_no_check(points : Vec<Point>) -> Self {
        if points.is_empty() { return Self { points : points, dim : 0 } }
        let dim = points[0].coords.len();
        Self { points : points, dim : dim }
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


    pub fn sum(&self, points : &Vec<usize>) -> Point {
        if points.len() == 0 {
            return Point::new(vec![0.0; self.dim()]);
        }

        let mut p = Point::new(Vec::new());
        for i in 0..self.dim() {
            let mut s = 0.0;
            for &x in points {
                s += self.get(x).coords[i];
            }
            p.coords.push(s);
        }
        p
    }

}
