


// bad practice having a utils file? might move / delete later
use crate::geometry::PointSet;
use crate::topology::Simplex;
use crate::construction::{ Adjacency, Distance };

use std::collections::HashMap;
use itertools::Itertools;


pub fn initialize_simplices(
    space : &PointSet,
    epsilon : f64,
    save_distance : bool
) -> (Vec<Simplex>, Adjacency, Distance) {
    let mut simplices : Vec<Simplex> = Vec::new();
    let mut adjacency : Adjacency = HashMap::new();
    let mut distance : Distance = HashMap::new();

    (0..space.len()).for_each(|x| simplices.push(Simplex::new(vec![x], 0.0))); // dim = 0

    for v in (0..space.len()).combinations(2) {
        let (x, y) = (v[0], v[1]); // x < y
        let d = space.get(x).distance(space.get(y));

        if d > epsilon { continue }
        simplices.push(Simplex::new(v, d)); // dim = 1

        if save_distance { distance.insert((x, y), d); }
        adjacency.entry(x).and_modify(|u| u.push(y)).or_insert(vec![y]);
        adjacency.entry(y).and_modify(|u| u.push(x)).or_insert(vec![x]);
    }

    //adjacency[i] already ordered for all i as per property of combinations
    (simplices, adjacency, distance)
}


pub fn sort_simplices(simplices : &mut Vec<Simplex>) {
    simplices.sort_by(|a, b| {
        a.filtration_value
            .partial_cmp(&b.filtration_value)
            .unwrap()
            .then(a.dim().cmp(&b.dim()))
            .then(a.vertices.cmp(&b.vertices))
    });
}
