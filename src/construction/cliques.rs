


use crate::geometry::PointSet;
use crate::construction::Construction;


// based on
// C. Bron, J. Kerbosch, Finding All Cliques of an Undirected Graph, 1973
// https://doi.org/10.1145/362342.362367
// NOTE similar to the algorithm described in https://ieeexplore.ieee.org/document/1559964
// TODO try : https://arxiv.org/abs/2311.13798v2
// assumes candidates and adjacency[i] are already sorted
pub fn bron_kerbosch(
    candidates : Vec<usize>,
    space : &PointSet,
    radius : fn(&Vec<usize>, &PointSet, &Construction) -> Option<f64>,
    cons : &mut Construction
) {
    // TODO add degeneracy ordering
    // do not forget that candidates still has to be ordered
    // as described in https://arxiv.org/abs/1006.5440
    // https://en.wikipedia.org/wiki/Degeneracy_(graph_theory)#Algorithms
    bron_kerbosch_rec(Vec::new(), candidates, space, radius, cons)
}


fn bron_kerbosch_rec(
    clique : Vec<usize>,
    candidates : Vec<usize>,
    space : &PointSet,
    radius : fn(&Vec<usize>, &PointSet, &Construction) -> Option<f64>,
    cons : &mut Construction
) {
    if clique.len() > 2 {
        match radius(&clique, space, cons) {
            Some(d) => cons.push(clique.clone(), d),
            None => return
        }
    }

    if clique.len() == cons.max_dim + 1 || candidates.len() == 0 { return; }
    if clique.len() + candidates.len() < 3 { return; }

    for (i, &x) in candidates.iter().enumerate() {
        bron_kerbosch_rec(
            join_back(&clique, x),
            intersection_ordered(&candidates[i..].to_vec(), &cons.adjacency[&x]),
            space, radius, cons
        )
    }
}


fn intersection_ordered(a : &Vec<usize>, b : &Vec<usize>) -> Vec<usize> {
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


fn join_back(v : &Vec<usize>, x : usize) -> Vec<usize> {
    [v, vec![x].as_slice()].concat()
}
