


use crate::geometry::PointSet;
use crate::topology::{ Simplex, Filtration };

use crate::construction::miniball;
use crate::construction::common::*;

use std::collections::HashMap;
use itertools::Itertools;


pub fn cech(point_set : &PointSet, max_epsilon : Option<f64>, max_dim : Option<usize>) -> Filtration {
    let max_epsilon = max_epsilon.unwrap_or(f64::MAX / 2.0);
    let max_dim = max_dim.unwrap_or(usize::MAX - 1);

    let mut simplices : Vec<Simplex> = Vec::new();
    (0..point_set.len()).for_each(|x| simplices.push(Simplex::new(vec![x], 0.0))); // dim = 0

    let mut adjacency : HashMap<usize, Vec<usize>> = HashMap::new();

    for v in (0..point_set.len()).combinations(2) {
        let (x, y) = (v[0], v[1]); // x < y
        let d = point_set.get(x).distance(point_set.get(y));

        if d > 2.0 * max_epsilon { continue }
        simplices.push(Simplex::new(v, d / 2.0)); // dim = 1

        adjacency.entry(x).and_modify(|u| u.push(y)).or_insert(vec![y]);
        adjacency.entry(y).and_modify(|u| u.push(x)).or_insert(vec![x]);
    }

    // TODO degeneracy ordering
    let candidates : Vec<usize> = (0..point_set.len()).collect(); // has to be ordered
    cliques(
        Vec::new(),
        candidates,
        max_dim + 1,
        &adjacency,
        point_set,
        max_epsilon,
        &mut simplices
    );

    simplices.sort_by(|a, b| {
        a.filtration_value
            .partial_cmp(&b.filtration_value)
            .unwrap()
            .then(a.dimension().cmp(&b.dimension()))
            .then(a.vertices.cmp(&b.vertices))
    });

    Filtration::new(simplices)
}


// same as for vietoris rips with some extra checks
fn cliques(
    clique : Vec<usize>,
    candidates : Vec<usize>,
    max_k : usize,
    adjacency : &HashMap<usize, Vec<usize>>,
    space : &PointSet,
    max_epsilon : f64,
    result : &mut Vec<Simplex>
) {
    if clique.len() > 2 {
        match in_ball(&clique, max_epsilon, space) {
            Some(d) => result.push(Simplex::new(clique.clone(), d)),
            None => return
        }
    }

    if clique.len() == max_k || candidates.len() == 0 { return; }
    if clique.len() + candidates.len() < 3 { return; }

    for (i, &x) in candidates.iter().enumerate() {
        cliques(
            [&clique, vec![x].as_slice()].concat(),
            intersection_ordered(&candidates[i..].to_vec(), &adjacency[&x]),
            max_k,
            adjacency,
            space,
            max_epsilon,
            result
        )
    }
}


//TODO change algorithm used based on dimesion and size of clique
fn in_ball(clique : &Vec<usize>, max_epsilon : f64, space : &PointSet) -> Option<f64> {
    let miniball = miniball::welzl(clique, space);
    if miniball.r() > max_epsilon { return None; }
    Some(miniball.r())
}
