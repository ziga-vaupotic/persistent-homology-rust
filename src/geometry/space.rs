use crate::geometry::Point;

pub trait Metric: Copy {
    fn distance(&self, a: &Point, b: &Point) -> f64;
}

pub trait Norm: Copy {
    fn norm(&self, a: &Point) -> f64;
}

pub trait InnerProduct: Copy {
    fn dot(&self, a: &Point, b: &Point) -> f64;
}

impl<T: Norm> Metric for T {
    fn distance(&self, a: &Point, b: &Point) -> f64 {
        let diff = a - b;
        self.norm(&diff)
    }
}

impl<T: InnerProduct> Norm for T {
    fn norm(&self, a: &Point) -> f64 {
        self.dot(a, a).sqrt()
    }
}

pub trait Euclidean: InnerProduct {} // specific inner product space

#[derive(Clone, Copy)]
pub struct EuclideanInnerProduct;

impl InnerProduct for EuclideanInnerProduct {
    fn dot(&self, a: &Point, b: &Point) -> f64 {
        a.coords.dot(&b.coords)
    }
}

impl Euclidean for EuclideanInnerProduct {}
