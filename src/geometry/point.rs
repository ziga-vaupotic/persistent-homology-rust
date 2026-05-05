#[derive(Debug, Clone, PartialEq)]
pub struct Point {
    pub coords: Vec<f64>,
}

impl Point {
    pub fn new(coords: Vec<f64>) -> Self {
        Self { coords }
    }

    pub fn dimension(&self) -> usize {
        self.coords.len()
    }

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
    let a = Point::new(vec![0.0, 0.0]);
    let b = Point::new(vec![3.0, 4.0]);
    assert_eq!(a.distance(&b), 5.0);
}