use crate::geometry::{Point};

pub trait Metric {
    fn distance(&self, a: &Point, b: &Point) -> f64;
}

pub trait InnerProduct {
    fn dot(&self, a: &Point, b: &Point) -> f64
}

impl<T: InnerProduct> Metric for T {
    fn distance(&self, a: &Point, b: &Point) -> f64 {
        let dif = a -b;
        self.dot(&diff, &diff).sqrt()
    }
}