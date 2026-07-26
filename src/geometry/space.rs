use crate::geometry::Point;

/// A metric space over points.
///
/// Implementations define a distance function satisfying the metric axioms.
/// The metric must be `Copy` so that point clouds can store the geometry by value.
pub trait Metric: Copy {
    /// Compute the distance between two points.
    fn distance(&self, a: &Point, b: &Point) -> f64;
}

/// A normed vector space over points.
///
/// Implements a norm (length) function on points.
pub trait Norm: Copy {
    /// Compute the norm (length) of a point/vector.
    fn norm(&self, a: &Point) -> f64;
}

/// An inner product space over points.
///
/// Implements a bilinear inner (dot) product function.
pub trait InnerProduct: Copy {
    /// Compute the inner (dot) product of two points/vectors.
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

/// Marker trait for Euclidean inner product spaces.
pub trait Euclidean: InnerProduct {}

#[derive(Clone, Copy)]
/// A Euclidean inner product implementation over points.
///
/// Computes the standard dot product: $\langle a, b \rangle = a \cdot b$.
/// Induces the Euclidean norm and metric.
pub struct EuclideanInnerProduct;

impl InnerProduct for EuclideanInnerProduct {
    fn dot(&self, a: &Point, b: &Point) -> f64 {
        a.coords.dot(&b.coords)
    }
}

impl Euclidean for EuclideanInnerProduct {}
