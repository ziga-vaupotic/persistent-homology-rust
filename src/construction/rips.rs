use crate::geometry::point::Point;
use crate::geometry::point_set::PointSet;

use crate::topology::simplex::Simplex;
use crate::topology::filtration::Filtration;

use itertools::Itertools;

pub fn vietoris_rips(points: &PointSet, max_dim: usize, max_epsilon: Option<f64>) -> Filtration {
    let max_epsilon = max_epsilon.unwrap_or(f64::MAX);

    let mut simplices = Vec::new();
    let n = points.len();


    // precompute distance in a distance matrix

    let mut dist = vec![vec![0.0; n]; n];

    for i in 0..n {
        for j in (i + 1)..n {
            let d = points.get(i).distance(&points.get(j));
            dist[i][j] = d;
            dist[j][i] = d;
        }
    }


    // Generate simplices of dimension 0 to max_dim
    for dim in 0..=max_dim {
        for combo in (0..n).combinations(dim + 1) {
        let mut max_dist: f64 = 0.0;

        let valid = combo.iter().enumerate().all(|(a, &i)| {
            combo[a + 1..].iter().all(|&j| {
                let d = dist[i][j];
                max_dist = max_dist.max(d);
                d <= max_epsilon
            })
        });

        if valid {
            simplices.push(Simplex::new(combo, max_dist));
        }
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

        let pointset = PointSet::new(points).expect("Pointset couldn't be generated.");

        let filtration = vietoris_rips(&pointset, 2, None);

        // Should have 3 0-simplices, 3 1-simplices, 1 2-simplex
        assert_eq!(filtration.simplices.len(), 7);

        // 0-simplices have filtration 0.0, others 0.5
        for i in 0..3 {
            assert_eq!(filtration.simplices[i].filtration_value, 0.0);
        }
        for i in 3..7 {
            assert!((filtration.simplices[i].filtration_value - 1.0).abs() < 1e-10);
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

        let pointset = PointSet::new(points).expect("Pointset couldn't be generated.");

        let filtration = vietoris_rips(&pointset, 2, None);

        assert_eq!(filtration.simplices.len(), 1);
        assert_eq!(filtration.simplices[0].vertices, vec![0]);
        assert_eq!(filtration.simplices[0].filtration_value, 0.0);
    }
}