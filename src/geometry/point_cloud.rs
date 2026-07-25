

use nalgebra::DVector;

use crate::geometry::{ Ball, Metric, InnerProduct };


#[derive(Clone)]
pub struct Point {
    pub coords: DVector<f64>,
}


impl Point {

    pub fn new(coords : Vec<f64>) -> Self {
        Self { coords : DVector::from_vec(coords) }
    }


    pub fn new_from_dvec(coords : DVector<f64>) -> Self {
        Self { coords : coords }
    }


    pub fn len(&self) -> usize {
        self.coords.len()
    }


    pub fn is_zero(&self) -> bool {
        self.coords.iter().all(|x| x.abs() <= 1e-14)
    }


    pub fn sum(points : &Vec<Point>) -> Result<Point, String> {
        if points.len() == 0 {
            return Err("no points given".into());
        }
        let dim = points[0].len();
        if !points.iter().all(|x| x.len() == dim) {
            return Err("inconsistent dimension".into());
        }

        let p = Point::new(Vec::new());
        for i in 0..dim {
            let mut s = 0.0;
            points.iter().for_each(|x| s += x.coords[i]);
            let _ = p.coords.push(s);
        }
        Ok(p)
    }


    pub fn sum_no_check(points : &Vec<Point>) -> Point {
        if points.len() == 0 {
            return Point::new(Vec::new());
        }

        let p = Point::new(Vec::new());
        for i in 0..points[0].len() {
            let mut s = 0.0;
            points.iter().for_each(|x| s += x.coords[i]);
            let _ = p.coords.push(s);
        }
        p
    }


    pub fn multiply(&mut self, lambda : f64) {
        self.coords = lambda * self.coords.clone();
    }


    pub fn standard_unit(i : usize, dim : usize) -> Self {
        let mut coords = vec![0.0; dim];
        coords[i] = 1.0;
        Point::new(coords)
    }

}


impl std::ops::Add<&Point> for &Point {
    type Output = Point;

    fn add(self, right : &Point) -> Point {
        let mut sum = Point::new(Vec::new());
        sum.coords = self.coords.clone() + right.coords.clone();
        sum
    }

}


impl std::ops::Sub<&Point> for &Point {
    type Output = Point;

    fn sub(self, right : &Point) -> Point {
        let mut sum = Point::new(Vec::new());
        sum.coords = self.coords.clone() - right.coords.clone();
        sum
    }

}


pub struct PointCloud <M> { // over R^n
    points: Vec<Point>,
    geometry : M, // Metric space e.g.
    dim: usize, // This is added as an enforcement of consistency in a nodeset
}


impl<M> PointCloud<M>
where M : Copy
{
    pub fn new(points: Vec<Point>, geometry: M) -> Result<Self, String> {
        if points.is_empty() {
            return Ok(Self { points, geometry : geometry, dim: 0 });
        }

        let dim = points[0].coords.len();
        if !points.iter().all(|p| p.coords.len() == dim) {
            return Err("Inconsistent point dimensions".into());
        }

        Ok(Self { points, geometry : geometry, dim : dim })
    }


    pub fn new_no_check(points : Vec<Point>, geometry : M) -> Self {
        if points.is_empty() { return Self { points, geometry : geometry, dim: 0 } }
        let dim = points[0].coords.len();
        Self { points, geometry : geometry, dim: dim }
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

    pub fn get_geometry(&self) -> M {
        self.geometry
    }

}


impl<M> PointCloud<M>
where 
    M: Metric,
{

    pub fn distance(&self, a : &Point, b : &Point) -> f64 {
        self.geometry.distance(a, b)
    }


    pub fn contained_in_ball(&self, ball: &Ball, point: &Point) -> bool {
        if ball.o().len() == 0 { return false; }
        self.distance(ball.o(), point) <= ball.r()
    }

}


impl<M> PointCloud<M>
where
    M: InnerProduct
{

    pub fn norm_squared(&self, a : &Point) -> f64 {
        self.geometry.dot(a, a)
    }


    pub fn norm(&self, a : &Point) -> f64 {
        self.norm_squared(a).sqrt()
    }

}
