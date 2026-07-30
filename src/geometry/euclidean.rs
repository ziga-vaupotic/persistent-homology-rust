use super::{InnerProductSpace, Space, VectorSpace};
use nalgebra::SVector;
use std::ops::{Add, Mul, Sub};

// TODO docs
#[derive(Clone)]
pub struct EuclideanSpace<const N: usize>;

impl<const N: usize> Space for EuclideanSpace<N> {
    type Element = Point<N>;
}

impl<const N: usize> VectorSpace for EuclideanSpace<N> {
    type Scalar = f64;

    fn add(a: &Self::Element, b: &Self::Element) -> Self::Element {
        Point::new(a.coords + b.coords)
    }

    fn sub(a: &Self::Element, b: &Self::Element) -> Self::Element {
        Point::new(a.coords - b.coords)
    }

    fn mul(a: &Self::Element, b: &Self::Scalar) -> Self::Element {
        Point::new(*b * a.coords)
    }

    fn dim() -> usize {
        N
    }
}

impl<const N: usize> InnerProductSpace for EuclideanSpace<N> {
    fn dot(a: &Self::Element, b: &Self::Element) -> f64 {
        a.coords.dot(&b.coords)
    }
}

impl<const N: usize> EuclideanSpace<N> {
    pub fn standard_unit(i: usize) -> Point<N> {
        assert!(i < N);
        let mut coords = SVector::zeros();
        coords[i] = 1.0;

        Point::new(coords)
    }

    pub fn zero() -> Point<N> {
        let coords = SVector::zeros();
        Point::new(coords)
    }
}

#[derive(Clone)]
pub struct Point<const N: usize> {
    pub coords: SVector<f64, N>,
}

impl<const N: usize> Point<N> {
    pub fn new(coords: impl Into<SVector<f64, N>>) -> Self {
        Self {
            coords: coords.into(),
        }
    }

    pub fn is_zero(&self) -> bool {
        self.coords.iter().all(|x| x.abs() <= 1e-14)
    }
}

impl<const N: usize> Add<&Point<N>> for &Point<N> {
    type Output = Point<N>;

    fn add(self, right: &Point<N>) -> Point<N> {
        EuclideanSpace::add(self, right)
    }
}

impl<const N: usize> Sub<&Point<N>> for &Point<N> {
    type Output = Point<N>;

    fn sub(self, right: &Point<N>) -> Point<N> {
        EuclideanSpace::sub(self, right)
    }
}

impl<const N: usize> Mul<&Point<N>> for f64 {
    type Output = Point<N>;

    fn mul(self, right: &Point<N>) -> Point<N> {
        EuclideanSpace::mul(right, &self)
    }
}
