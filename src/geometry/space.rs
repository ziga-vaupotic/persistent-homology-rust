// TODO docs
pub trait Set {
    type Element: Clone;
}

pub trait VectorSpace: Set {
    type Scalar;

    fn add(a: &Self::Element, b: &Self::Element) -> Self::Element;
    fn sub(a: &Self::Element, b: &Self::Element) -> Self::Element;
    fn mul(a: &Self::Element, b: &Self::Scalar) -> Self::Element;
    fn dim() -> usize;
}

pub trait Space {
    type Set: Set;
}

pub trait MetricSpace: Space {
    fn distance(a: &<Self::Set as Set>::Element, b: &<Self::Set as Set>::Element) -> f64;
}

pub trait NormedSpace: Space<Set: VectorSpace> {
    fn norm(a: &<Self::Set as Set>::Element) -> f64;
}

impl<T: NormedSpace> MetricSpace for T {
    fn distance(a: &<Self::Set as Set>::Element, b: &<Self::Set as Set>::Element) -> f64 {
        Self::norm(&Self::Set::sub(a, b))
    }
}

pub trait InnerProductSpace: Space<Set: VectorSpace> {
    fn dot(a: &<Self::Set as Set>::Element, b: &<Self::Set as Set>::Element) -> f64;
}

impl<T: InnerProductSpace> NormedSpace for T {
    fn norm(a: &<Self::Set as Set>::Element) -> f64 {
        Self::dot(a, a).sqrt()
    }
}
