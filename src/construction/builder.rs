use crate::geometry::{Metric, PointCloud};
use crate::topology::Simplex;

use itertools::Itertools;
use std::collections::HashMap;

pub struct Construction {
    pub simplices: Vec<Simplex>,
    pub adjacency: HashMap<usize, Vec<usize>>, // adjacency[v] = N(v)
    pub distance: HashMap<(usize, usize), f64>,

    pub max_k: usize,
    pub max_epsilon: f64,
    pub tolerance: f64,
}

impl Construction {
    pub fn new<M>(max_dim: usize, max_epsilon: f64, tolerance: f64, space: &PointCloud<M>) -> Self
    where
        M: Metric,
    {
        let mut simplices: Vec<Simplex> = Vec::new();
        let mut adjacency: HashMap<usize, Vec<usize>> = HashMap::new();
        let distance: HashMap<(usize, usize), f64> = HashMap::new();

        let n = space.len();
        (0..n).for_each(|x| simplices.push(Simplex::new(vec![x], 0.0))); // dim = 0
        (0..n).for_each(|x| {
            adjacency.insert(x, Vec::new());
        });

        Self {
            simplices,
            adjacency,
            distance,
            max_k: max_dim + 1,
            max_epsilon,
            tolerance,
        }
    }

    pub fn push(&mut self, clique: Vec<usize>, filtration_value: f64) {
        self.simplices.push(Simplex::new(clique, filtration_value));
    }

    pub fn traverse_edges<M>(
        &mut self,
        space: &PointCloud<M>,
        factor: f64,
        save_distance: bool,
    ) -> bool
    where
        M: Metric,
    {
        let n = space.len();

        let mut has_edges = false;
        for v in (0..n).combinations(2) {
            let (x, y) = (v[0], v[1]); // x < y
            let d = space.distance(space.get(x), space.get(y));

            if d > factor * self.max_epsilon {
                continue;
            }
            self.push(v, d / factor); // dim = 1

            if save_distance {
                self.distance.insert((x, y), d);
            }
            self.adjacency.entry(x).and_modify(|u| u.push(y));
            self.adjacency.entry(y).and_modify(|u| u.push(x));

            has_edges = true;
        }
        //adjacency[i] already ordered for all i as per property of combinations
        has_edges
    }

    pub fn sort_simplices(&mut self) {
        self.simplices.sort_by(|a, b| {
            a.filtration_value
                .partial_cmp(&b.filtration_value)
                .unwrap()
                .then(a.dim().cmp(&b.dim()))
                .then(a.vertices.cmp(&b.vertices))
        });
    }
}
