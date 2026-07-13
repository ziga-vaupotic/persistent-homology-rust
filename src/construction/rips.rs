


use crate::geometry::PointSet;

use crate::topology::{ Simplex, Filtration };

use std::collections::HashMap;
use itertools::Itertools;

use crate::construction::common::*;


/*
// for graph density
fn binom(a : usize, b : usize) -> usize {
    if b == 0 { return 1; }
    (a * binom(a - 1, b - 1)) / b
}
*/

fn rips_simplex( clique : Vec<usize>, distance : &HashMap<(usize, usize), f64>) -> Simplex {
    let mut max_d = 0.0;
    for mut v in (0..clique.len()).combinations(2) {
        v.sort();
        let &d = distance.get(&(clique[v[0]], clique[v[1]])).unwrap();
        max_d = if d >= max_d { d } else { max_d };
    }
    Simplex::new(clique, max_d)
}


// using modified Bron-Kerbosch (1973)
// similar to as described in : https://ieeexplore.ieee.org/document/1559964
// TODO try : https://arxiv.org/abs/2311.13798v2
fn cliques(
    clique : Vec<usize>,
    candidates : Vec<usize>,
    max_k : usize,
    adjacency : &HashMap<usize, Vec<usize>>,
    distance : &HashMap<(usize, usize), f64>,
    result : &mut Vec<Simplex>
) {
    if clique.len() >= 3 { result.push(rips_simplex(clique.clone(), distance)); }

    if clique.len() == max_k || candidates.len() == 0 { return; }
    if clique.len() + candidates.len() < 3 { return; }

    for (i, &x) in candidates.iter().enumerate() {
        cliques(
            [&clique, vec![x].as_slice()].concat(),
            intersection_ordered(&candidates[i..].to_vec(), &adjacency[&x]),
            max_k,
            adjacency,
            distance,
            result
        )
    }
}


pub fn vietoris_rips(point_set : &PointSet, max_epsilon : Option<f64>, max_dim : Option<usize>) -> Filtration {
    let max_epsilon = max_epsilon.unwrap_or(f64::MAX);
    let max_dim = max_dim.unwrap_or(usize::MAX - 1);

    let mut simplices : Vec<Simplex> = Vec::new();
    (0..point_set.len()).for_each(|x| simplices.push(Simplex::new(vec![x], 0.0))); // dim = 0

    //let mut num_edges = 0;
    let mut adjacency : HashMap<usize, Vec<usize>> = HashMap::new();
    let mut distance : HashMap<(usize, usize), f64> = HashMap::new();

    for v in (0..point_set.len()).combinations(2) {
        let (x, y) = (v[0], v[1]); // x < y
        let d = point_set.get(x).distance(point_set.get(y));

        if d > max_epsilon { continue }
        simplices.push(Simplex::new(v, d)); // dim = 1
        //num_edges += 1;

        distance.insert((x, y), d);
        adjacency.entry(x).and_modify(|u| u.push(y)).or_insert(vec![y]);
        adjacency.entry(y).and_modify(|u| u.push(x)).or_insert(vec![x]);
    }
    //adjacency[i] already ordered for all i as per property of combinations

    if adjacency.is_empty() { return Filtration::new(simplices) }

    /*
    let len_choose_2 = binom(point_set.len(), 2);
    let f : f64 = num_edges as f64 / len_choose_2 as f64;
    println!("graph density {} / {} = {:?}", num_edges, len_choose_2, f);
    */

    // TODO add degeneracy ordering
    // order by degeneracy get some permutation of 0..n,
    // input candidates = 0..n with phi : 0..n -> 0..n bijection that maps i to element at index i
    // in degeneracy ordering
    // as described in : https://arxiv.org/abs/1006.5440
    let candidates : Vec<usize> = (0..point_set.len()).collect(); // has to be ordered
    cliques(
        Vec::new(),
        candidates,
        max_dim + 1,
        &adjacency,
        &distance,
        &mut simplices
    );

    simplices.sort_by(|a, b| {
        a.filtration_value
            .partial_cmp(&b.filtration_value)
            .unwrap()
            .then(a.dim().cmp(&b.dim()))
            .then(a.vertices.cmp(&b.vertices))
    });

    Filtration::new(simplices)
}


#[cfg(test)]
mod tests {

    use super::*;


    #[test]
    fn _test_vietoris_rips_three_points() {
        use std::f64::consts::PI;
        let points = vec![
            Point::new(vec![f64::cos(0.0), f64::sin(0.0)]),
            Point::new(vec![f64::cos(2.0 * PI / 3.0), f64::sin(2.0 * PI / 3.0)]),
            Point::new(vec![f64::cos(4.0 * PI / 3.0), f64::sin(4.0 * PI / 3.0)])
        ];
        let pointset = PointSet::new(points).expect("Pointset couldn't be generated.");

        let filtration = vietoris_rips(&pointset, None, None);

        // should have [0], [1], [2], [0, 1], [0, 2], [1, 2], [0, 1, 2]
        assert_eq!(filtration.simplices.len(), 7);

        // 0-simplices have filtration 0.0, 1 and 2-simplices 2 cos(PI / 6)
        for i in 0..3 {
            assert_eq!(filtration.simplices[i].filtration_value, 0.0);
        }
        for i in 3..7 {
            assert!((filtration.simplices[i].filtration_value - 2.0 * f64::cos(PI / 6.0)).abs() < 1e-12);
        }

        let dims: Vec<usize> = filtration.simplices.iter().map(|s| s.dim()).collect();
        assert_eq!(dims, vec![0, 0, 0, 1, 1, 1, 2]);

        assert_eq!(filtration.simplices[6].vertices, vec![0, 1, 2]);
    }

    #[test]
    fn test_vietoris_rips_single_point() {
        let points = vec![Point::new(vec![0.0, 0.0])];
        let pointset = PointSet::new(points).expect("Pointset couldn't be generated.");

        let filtration = vietoris_rips(&pointset, None, None);

        assert_eq!(filtration.simplices.len(), 1);
        assert_eq!(filtration.simplices[0].vertices, vec![0]);
        assert_eq!(filtration.simplices[0].filtration_value, 0.0);
    }
}
