// TODO docs
pub trait Space {
    type Element;
}

pub trait MetricSpace: Space {
    fn distance(a: &Self::Element, b: &Self::Element) -> f64;
}

pub trait VectorSpace: Space {
    type Scalar;

    fn add(a: &Self::Element, b: &Self::Element) -> Self::Element;
    fn sub(a: &Self::Element, b: &Self::Element) -> Self::Element;
    fn mul(a: &Self::Element, b: &Self::Scalar) -> Self::Element;
    fn dim() -> usize;
}

pub trait NormSpace: VectorSpace {
    fn norm(a: &Self::Element) -> f64;
}

impl<T: NormSpace> MetricSpace for T {
    fn distance(a: &Self::Element, b: &Self::Element) -> f64 {
        Self::norm(&Self::sub(a, b))
    }
}

pub trait InnerProductSpace: VectorSpace<Scalar = f64> {
    fn dot(a: &Self::Element, b: &Self::Element) -> Self::Scalar;
    fn norm_squared(a: &Self::Element) -> f64 {
        Self::dot(a, a)
    }
}

impl<T: InnerProductSpace> NormSpace for T {
    fn norm(a: &Self::Element) -> f64 {
        Self::norm_squared(a).sqrt()
    }
}
