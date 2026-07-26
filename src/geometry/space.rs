use crate::geometry::Point;

/// A metric space over points.
///
/// Implementations define a distance function satisfying the metric axioms.
/// The metric must be `Copy` so that point clouds can store the geometry by value.
pub trait Metric<const D: usize>: Copy {
    /// Compute the distance between two points.
    fn distance(&self, a: &Point<D>, b: &Point<D>) -> f64;
}

/// A normed vector space over points.
///
/// Implements a norm (length) function on points.
pub trait Norm<const D: usize>: Copy {
    /// Compute the norm (length) of a point/vector.
    fn norm(&self, a: &Point<D>) -> f64;
}

/// An inner product space over points.
///
/// Implements a bilinear inner (dot) product function.
pub trait InnerProduct<const D: usize>: Copy {
    /// Compute the inner (dot) product of two points/vectors.
    fn dot(&self, a: &Point<D>, b: &Point<D>) -> f64;
}

impl<T, const D: usize> Metric<D> for T
where
    T: Norm<D>,
{
    fn distance(&self, a: &Point<D>, b: &Point<D>) -> f64 {
        let diff = a - b;
        self.norm(&diff)
    }
}

impl<T, const D: usize> Norm<D> for T
where
    T: InnerProduct<D>,
{
    fn norm(&self, a: &Point<D>) -> f64 {
        self.dot(a, a).sqrt()
    }
}

/// Marker trait for Euclidean inner product spaces.
pub trait Euclidean<const D: usize>: InnerProduct<D> {}

#[derive(Clone, Copy)]
/// A Euclidean inner product implementation over points.
///
/// Computes the standard dot product: $\langle a, b \rangle = a \cdot b$.
/// Induces the Euclidean norm and metric.
pub struct EuclideanInnerProduct;

impl<const D: usize> InnerProduct<D> for EuclideanInnerProduct {
    fn dot(&self, a: &Point<D>, b: &Point<D>) -> f64 {
        a.coords.dot(&b.coords)
    }
}

impl<const D: usize> Euclidean<D> for EuclideanInnerProduct {}
