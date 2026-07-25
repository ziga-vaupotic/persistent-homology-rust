


use crate::geometry::space::*;
use crate::geometry::{Ball};


#[derive(Clone)]
pub struct Point {
    pub coords: Vec<f64>,
}


impl Point {
    pub fn new(coords : Vec<f64>) -> Self {
        Self { coords : coords }
    }


    pub fn len(&self) -> usize {
        self.coords.len()
    }
}


pub struct PointCloud <M> {
    points: Vec<Point>,
    metric: M, // Metric space e.g.
    dim: usize, // This is added as an enforcement of consistency in a nodeset
}


impl<M> PointCloud<M>
{
    pub fn new(points: Vec<Point>, metric: M) -> Result<Self, String> {
        if points.is_empty() {
            return Ok(Self { points, metric, dim: 0 });
        }

        let dim = points[0].coords.len();
        if !points.iter().all(|p| p.coords.len() == dim) {
            return Err("Inconsistent point dimensions".into());
        }

        Ok(Self { points, metric, dim })
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
            points.iter().for_each(|&x| s += self.get(x).coords[i]);
            p.coords.push(s);
        }
        p
    }
}


impl<M> PointCloud<M>
where 
    M: Metric,
{

    pub fn contained_in_ball(&self, ball: &Ball, point: &Point) -> bool {
        self.metric.distance(ball.o(), point) <= ball.r()
    }

}
