use crate::geometry::point::Point;
use crate::topology::simplex::Simplex;
use crate::topology::filtration::Filtration;

use itertools::Itertools;

pub fn vietoris_rips(points: &[Point], max_dim: usize) -> Filtration {
    let mut simplices = Vec::new();
    let n = points.len();

    // Generate simplices of dimension 0 to max_dim
    for dim in 0..=max_dim {
        for combo in (0..n).combinations(dim + 1) {
            // Compute max pairwise distance
            let mut max_dist = 0.0;
            for &i in &combo {
                for &j in &combo {
                    if i < j {
                        let d = points[i].distance(&points[j]);
                        if d > max_dist {
                            max_dist = d;
                        }
                    }
                }
            }
            let filtration_value = max_dist / 2.0;
            let simplex = Simplex::new(combo, filtration_value);
            simplices.push(simplex);
        }
    }

    // Sort simplices by filtration_value, then by dimension, then by vertices
    simplices.sort_by(|a, b| {
        a.filtration_value
            .partial_cmp(&b.filtration_value)
            .unwrap()
            .then(a.dimension().cmp(&b.dimension()))
            .then(a.vertices.cmp(&b.vertices))
    });

    Filtration { simplices }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_vietoris_rips_three_points() {
        let points = vec![
            Point::new(vec![0.0, 0.0]),
            Point::new(vec![1.0, 0.0]),
            Point::new(vec![0.5, (3.0_f64).sqrt() / 2.0]),
        ];

        let filtration = vietoris_rips(&points, 2);

        // Should have 3 0-simplices, 3 1-simplices, 1 2-simplex
        assert_eq!(filtration.simplices.len(), 7);

        // 0-simplices have filtration 0.0, others 0.5
        for i in 0..3 {
            assert_eq!(filtration.simplices[i].filtration_value, 0.0);
        }
        for i in 3..7 {
            assert!((filtration.simplices[i].filtration_value - 0.5).abs() < 1e-10);
        }

        // Check dimensions
        let dims: Vec<usize> = filtration.simplices.iter().map(|s| s.dimension()).collect();
        assert_eq!(dims, vec![0, 0, 0, 1, 1, 1, 2]);

        // Check vertices for 2-simplex
        assert_eq!(filtration.simplices[6].vertices, vec![0, 1, 2]);
    }

    #[test]
    fn test_vietoris_rips_single_point() {
        let points = vec![Point::new(vec![0.0, 0.0])];
        let filtration = vietoris_rips(&points, 2);

        assert_eq!(filtration.simplices.len(), 1);
        assert_eq!(filtration.simplices[0].vertices, vec![0]);
        assert_eq!(filtration.simplices[0].filtration_value, 0.0);
    }
}