use crate::geometry::{Metric, PointCloud};
use crate::topology::Simplex;

use itertools::Itertools;
use std::collections::HashMap;

/// Intermediate state for constructing simplicial complexes.
///
/// The `Construction` struct accumulates simplices and maintains adjacency information
/// needed by the clique enumeration algorithm. It tracks distances and adjacency for efficient
/// complex construction.
pub struct Construction {
    /// Simplices accumulated during construction.
    pub simplices: Vec<Simplex>,
    /// Adjacency list: `adjacency[v]` contains all vertices adjacent to vertex `v`.
    pub adjacency: HashMap<usize, Vec<usize>>,
    /// Cache of pairwise distances between vertices (if requested).
    pub distance: HashMap<(usize, usize), f64>,

    /// Maximum simplex dimension (1 + number of vertices in largest simplex).
    pub max_k: usize,
    /// Maximum allowed filtration value.
    pub max_epsilon: f64,
    /// Tolerance for radius computations (used in Čech complex approximations).
    pub tolerance: f64,
}

impl Construction {
    /// Initialize a new construction from a point cloud.
    ///
    /// Creates vertex simplices (0-simplices) for each point and initializes the adjacency structure.
    ///
    /// # Arguments
    ///
    /// * `max_dim` - Maximum simplex dimension (number of vertices - 1).
    /// * `max_epsilon` - Maximum filtration value threshold.
    /// * `tolerance` - Tolerance for geometric computations.
    /// * `space` - The point cloud.
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

    /// Add a simplex to the construction.
    ///
    /// # Arguments
    ///
    /// * `clique` - Vertex indices of the simplex.
    /// * `filtration_value` - Filtration value at which this simplex appears.
    pub fn push(&mut self, clique: Vec<usize>, filtration_value: f64) {
        self.simplices.push(Simplex::new(clique, filtration_value));
    }

    /// Traverse and add all edges within the filtration threshold.
    ///
    /// Finds all pairs of vertices with distance at most `max_epsilon * factor`
    /// and adds them as 1-simplices. Optionally caches pairwise distances.
    ///
    /// # Arguments
    ///
    /// * `space` - The point cloud.
    /// * `factor` - Scaling factor for the distance threshold.
    /// * `save_distance` - Whether to cache pairwise distances.
    ///
    /// # Returns
    ///
    /// `true` if at least one edge was found, `false` if the graph is empty.
    ///
    /// # Note
    ///
    /// This method maintains the adjacency lists in sorted order for compatibility
    /// with the clique enumeration algorithm.
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

    /// Sort simplices by filtration value, dimension, and vertex indices.
    ///
    /// Arranges simplices in the order required for the persistence algorithm:
    /// first by increasing filtration value, then by dimension, and finally by vertex indices.
    /// This establishes a total order on simplices suitable for homology computation.
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
