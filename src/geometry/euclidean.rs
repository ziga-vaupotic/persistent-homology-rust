


use super::{ InnerProductSpace, Space, VectorSpace };
use std::ops::{ Add, Sub, Mul };


// TODO docs
pub struct EuclideanSpace<const N : usize>;


impl<const N : usize> Space for EuclideanSpace<N> {

    type Element = Point<N>;

}


impl<const N : usize> VectorSpace for EuclideanSpace<N> {
    type Scalar = f64;

    fn add(a : &Self::Element, b : &Self::Element) -> Self::Element {
        let mut coords = [0.0; N];
        (0..N).for_each(|i| coords[i] = a.coords[i] + b.coords[i]);
        Point { coords }
    }


    fn sub(a : &Self::Element, b : &Self::Element) -> Self::Element {
        let mut coords = [0.0; N];
        (0..N).for_each(|i| coords[i] = a.coords[i] - b.coords[i]);
        Point { coords }
    }


    fn mul(a : &Self::Element, b : &Self::Scalar) -> Self::Element {
        let mut coords = [0.0; N];
        (0..N).for_each(|i| coords[i] = a.coords[i] * b);
        Point { coords }
    }


    fn dim() -> usize {
        N
    }

}


impl<const N : usize> InnerProductSpace for EuclideanSpace<N> {

    fn dot(a : &Self::Element, b : &Self::Element) -> f64 {
        a.coords.iter().zip(b.coords.iter()).map(|(x, y)| x * y).sum()
    }

}


pub struct Point<const N : usize> {
    pub coords : [f64; N]
}


impl<const N : usize> Point<N> {

    pub fn from_fn(phi : fn(usize) -> f64) -> Self {
        let mut coords = [0.0; N];
        (0..N).for_each(|i| coords[i] = phi(i));
        Self { coords }
    }


    pub fn standard_unit(i : usize) -> Self {
        assert!(i < N);
        let mut coords = [0.0; N];
        coords[i] = 1.0;
        Self { coords }
    }


    pub fn is_zero(&self) -> bool {
        self.coords.iter().all(|x| x.abs() <= 1e-14)
    }

}


impl<const N : usize> Add<&Point<N>> for &Point<N> {
    type Output = Point<N>;

    fn add(self, right : &Point<N>) -> Point<N> {
        EuclideanSpace::add(self, right)
    }
}


impl<const N : usize> Sub<&Point<N>> for &Point<N> {
    type Output = Point<N>;

    fn sub(self, right : &Point<N>) -> Point<N> {
        EuclideanSpace::sub(self, right)
    }
}


impl<const N : usize> Mul<&Point<N>> for f64 {
    type Output = Point<N>;

    fn mul(self, right : &Point<N>) -> Point<N> {
        EuclideanSpace::mul(right, &self)
    }
}
