


#[derive(Debug, Clone, PartialEq)]
pub struct Point {
    pub coords: Vec<f64>,
}


impl Point {

    pub fn new(coords : Vec<f64>) -> Self {
        Self { coords : coords }
    }


    pub fn dim(&self) -> usize {
        self.coords.len()
    }


    pub fn is_zero(&self) -> bool {
        self.coords.iter().all(|x| x.abs() < 1e-12) // tolerance
    }


    pub fn dot(&self, other : &Self) -> f64 {
        self.coords.iter().zip(other.coords.iter()).map(|(x, y)| x * y).sum()
    }


    pub fn norm_square(&self) -> f64 {
        self.dot(self)
    }


    pub fn norm(&self) -> f64 {
        self.norm_square().sqrt()
    }


    pub fn normalize(&mut self) { // opinions on normalize / normalise conundrum???
        if self.is_zero() { return; }
        self.multiply(1.0 / self.norm());
    }


    pub fn distance_square(&self, other : &Self) -> f64 {
        self.coords.iter().zip(other.coords.iter()).map(|(x, y)| (x - y).powi(2)).sum::<f64>()
    }


    pub fn distance(&self, other : &Self) -> f64 {
        self.distance_square(other).sqrt()
    }


    pub fn multiply(&mut self, lambda : f64) {
        self.coords = self.coords.iter().map(|&x| lambda * x).collect();
    }


    pub fn add(&mut self, other : &Self) {
        for i in 0..self.dim() {
            self.coords[i] += other.coords[i];
        }
    }


    pub fn sum(points : &Vec<Point>) -> Result<Point, String> {
        if points.len() == 0 {
            return Err("no points given".into());
        }
        let dim = points[0].dim();
        if !points.iter().all(|x| x.dim() == dim) {
            return Err("inconsistent dimension".into());
        }

        let mut p = Point::new(Vec::new());
        for i in 0..dim {
            let mut s = 0.0;
            for x in points {
                s += x.coords[i];
            }
            p.coords.push(s);
        }
        Ok(p)
    }


    pub fn sum_no_check(points : &Vec<Point>) -> Point {
        if points.len() == 0 {
            return Point::new(Vec::new());
        }

        let mut p = Point::new(Vec::new());
        for i in 0..points[0].dim() {
            let mut s = 0.0;
            for x in points {
                s += x.coords[i];
            }
            p.coords.push(s);
        }
        p
    }


    pub fn subtract(&mut self, other : &Self) {
        for i in 0..self.dim() {
            self.coords[i] -= other.coords[i];
        }
    }


    // TODO add a safe version and move this to difference_no_check
    pub fn difference(a : &Point, b : &Point) -> Point {
        let mut dif : Vec<f64> = Vec::new();
        for i in 0..a.dim() {
            dif.push(a.coords[i] - b.coords[i]);
        }
        Point::new(dif)
    }


    pub fn project(&mut self, other : &Self) {
        let mut projection = other.clone();
        projection.multiply(self.dot(other) / other.norm_square());
        self.coords = projection.coords;
    }


    pub fn project_normal(&mut self, other : &Self) { // assumes other is normal
        let mut projection = other.clone();
        projection.multiply(self.dot(other));
        self.coords = projection.coords;
    }


    pub fn projection(u : &Point, v : &Point) -> Point { // returns a projection of u onto v
        let mut projection = u.clone();
        projection.project(v);
        projection
    }


    pub fn projection_normal(u : &Point, v : &Point) -> Point {
        let mut projection = u.clone();
        projection.project_normal(v);
        projection
    }


    pub fn standard_unit(i : usize, dim : usize) -> Point {
        let mut coords = vec![0.0; dim];
        coords[i] = 1.0;
        Point::new(coords)
    }

}


#[test]
fn test_distance_2d() {
    let a = Point::new(vec![0.0, 0.0]);
    let b = Point::new(vec![3.0, 4.0]);
    assert_eq!(a.distance(&b), 5.0);
}
