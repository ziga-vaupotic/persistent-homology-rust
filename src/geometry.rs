#[derive(Debug, Clone, PartialEq)]
pub struct Point<const D: usize> {
    pub coords: [f64; D],
}

impl<const D: usize> Point<D> {
    pub fn distance(&self, other: &Self) -> f64 {
        self.coords
            .iter()
            .zip(other.coords.iter())
            .map(|(a, b)| (a - b).powi(2))
            .sum::<f64>()
            .sqrt()
    }
}


#[test]
fn test_distance_2d() {
    let a = Point::<2> { coords: [0.0, 0.0] };
    let b = Point::<2> { coords: [3.0, 4.0] };
    assert_eq!(a.distance(&b), 5.0);
}