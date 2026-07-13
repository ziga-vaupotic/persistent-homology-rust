


use crate::geometry::PointSet;
use crate::topology::Simplex;
use crate::construction::{ Adjacency, Distance };


// based on
// C. Bron, J. Kerbosch, Finding All Cliques of an Undirected Graph, 1973
// https://doi.org/10.1145/362342.362367
// NOTE similar to the algorithm described in https://ieeexplore.ieee.org/document/1559964
// TODO try : https://arxiv.org/abs/2311.13798v2
pub fn bron_kerbosch(
    clique : Vec<usize>, // accepting recommendations to fix this mess of arguments
    candidates : Vec<usize>,
    max_k : usize,
    max_epsilon : f64,
    tolerance : f64,
    space : &PointSet,
    adjacency : &Adjacency,
    distance : &Distance,
    radius : fn(&Vec<usize>, f64, f64, &Distance, &PointSet) -> Option<f64>, // especially this
    result : &mut Vec<Simplex>
) {
    if clique.len() > 2 {
        match radius(&clique, max_epsilon, tolerance, distance, space) {
            Some(d) => result.push(Simplex::new(clique.clone(), d)),
            None => return
        }
    }

    if clique.len() == max_k || candidates.len() == 0 { return; }
    if clique.len() + candidates.len() < 3 { return; }

    for (i, &x) in candidates.iter().enumerate() {
        bron_kerbosch(
            [&clique, vec![x].as_slice()].concat(),
            intersection_ordered(&candidates[i..].to_vec(), &adjacency[&x]),
            max_k,
            max_epsilon,
            tolerance,
            space,
            adjacency,
            distance,
            radius,
            result
        )
    }
}


pub fn intersection_ordered(a : &Vec<usize>, b : &Vec<usize>) -> Vec<usize> {
    let (m, n) = (a.len(), b.len());
    let (mut i, mut j) = (0, 0);
    let mut intersection : Vec<usize> = Vec::new();
    while i < m && j < n {
        if a[i] < b[j] {
            i += 1;
            continue
        }
        if b[j] < a[i] {
            j += 1;
            continue
        }
        intersection.push(a[i]);
        i += 1;
        j += 1;
    }
    intersection
}
