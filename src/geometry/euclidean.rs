use super::{Ball, InnerProductSpace, PointCloud, Set, Space, VectorSpace};
use nalgebra::SVector;
use std::ops::{Add, Mul, Sub};

pub type EuclideanBall<const N: usize> = Ball<EuclideanSpace<N>>;
pub type EuclideanCloud<const N: usize> = PointCloud<EuclideanSpace<N>>;

// TODO docs
#[derive(Clone, Copy)]
pub struct EuclideanSpace<const N: usize>;

impl<const N: usize> Space for EuclideanSpace<N> {
    type Set = RealCoordinateSpace<N>;
}

impl<const N: usize> InnerProductSpace for EuclideanSpace<N> {
    fn dot(a: &Point<N>, b: &Point<N>) -> f64 {
        a.coords.dot(&b.coords)
    }
}

#[derive(Clone, Copy)]
pub struct RealCoordinateSpace<const N: usize>;

impl<const N: usize> Set for RealCoordinateSpace<N> {
    type Element = Point<N>;
}

impl<const N: usize> VectorSpace for RealCoordinateSpace<N> {
    type Scalar = f64;

    fn add(a: &Point<N>, b: &Point<N>) -> Point<N> {
        a + b
    }

    fn sub(a: &Point<N>, b: &Point<N>) -> Point<N> {
        a - b
    }

    fn mul(a: &Point<N>, b: &f64) -> Point<N> {
        *b * a
    }

    fn dim() -> usize {
        N
    }
}

impl<const N: usize> RealCoordinateSpace<N> {
    pub fn standard_unit(i: usize) -> Point<N> {
        assert!(i < N);
        let mut coords = SVector::zeros();
        coords[i] = 1.0;

        Point::new(coords)
    }

    pub fn zero() -> Point<N> {
        Point::new(SVector::zeros())
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
        Point::new(self.coords + right.coords)
    }
}

impl<const N: usize> Sub<&Point<N>> for &Point<N> {
    type Output = Point<N>;

    fn sub(self, right: &Point<N>) -> Point<N> {
        Point::new(self.coords - right.coords)
    }
}

impl<const N: usize> Mul<&Point<N>> for f64 {
    type Output = Point<N>;

    fn mul(self, right: &Point<N>) -> Point<N> {
        Point::new(self * right.coords)
    }
}
