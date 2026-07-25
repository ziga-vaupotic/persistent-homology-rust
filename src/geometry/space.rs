use crate::geometry::{Point};

pub trait Metric : Copy {
    fn distance(&self, a: &Point, b: &Point) -> f64;
}

pub trait InnerProduct : Copy {
    fn dot(&self, a: &Point, b: &Point) -> f64;
}

impl<T: InnerProduct> Metric for T {
    fn distance(&self, a: &Point, b: &Point) -> f64 {
        let diff = a - b;
        self.dot(&diff, &diff).sqrt()
    }
}


pub trait Euclidian : InnerProduct {}

#[derive(Clone, Copy)]
pub struct EuclidianInnerProduct;


impl InnerProduct for EuclidianInnerProduct {

    fn dot(&self, a : &Point, b : &Point) -> f64 {
        a.coords.dot(&b.coords)
    }

}


impl Euclidian for EuclidianInnerProduct {}
