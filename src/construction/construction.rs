


use crate::geometry::PointSet;
use crate::topology::Simplex;

use std::collections::HashMap;
use itertools::Itertools;


pub struct Construction {
    pub simplices : Vec<Simplex>,
    pub adjacency : HashMap<usize, Vec<usize>>,
    pub distance : HashMap<(usize, usize), f64>,

    pub max_dim : usize,
    pub max_epsilon : f64,
    pub tolerance : f64
}


impl Construction {

    pub fn new(max_dim : usize, max_epsilon : f64, tolerance : f64) -> Self {
        let simplices : Vec<Simplex> = Vec::new();
        let adjacency : HashMap<usize, Vec<usize>> = HashMap::new();
        let distance : HashMap<(usize, usize), f64> = HashMap::new();

        Self {
            simplices : simplices, adjacency : adjacency, distance : distance,
            max_dim : max_dim, max_epsilon : max_epsilon, tolerance : tolerance
        }
    }


    pub fn push(&mut self, clique : Vec<usize>, filtration_value : f64) {
        self.simplices.push(Simplex::new(clique, filtration_value));
    }


    pub fn traverse_edges(&mut self, space : &PointSet, epsilon : f64, save_distance : bool) {
        (0..space.len()).for_each(|x| self.push(vec![x], 0.0)); // dim = 0

        for v in (0..space.len()).combinations(2) {
            let (x, y) = (v[0], v[1]); // x < y
            let d = space.get(x).distance(space.get(y));

            if d > epsilon { continue }
            self.push(v, d); // dim = 1

            if save_distance { self.distance.insert((x, y), d); }
            self.adjacency.entry(x).and_modify(|u| u.push(y)).or_insert(vec![y]);
            self.adjacency.entry(y).and_modify(|u| u.push(x)).or_insert(vec![x]);
        }
        //adjacency[i] already ordered for all i as per property of combinations
    }

    pub fn no_adjacency(&self) -> bool {
        self.adjacency.is_empty()
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
