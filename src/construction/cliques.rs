use crate::construction::Construction;
use crate::geometry::{Metric, PointCloud};

pub fn find_all<M>(
    space: &PointCloud<M>,
    radius: fn(&[usize], &PointCloud<M>, &Construction) -> Option<f64>,
    cons: &mut Construction,
) where
    M: Metric,
{
    // TODO add degeneracy ordering
    // do not forget that candidates still has to be ordered
    // https://arxiv.org/abs/1006.5440
    // algorithm described in the article runs O(n d 3^(d / 3)) where d is degeneracy
    // https://en.wikipedia.org/wiki/Degeneracy_(graph_theory)#Algorithms
    // NOTE somewhat unviable with the current setup --- would need major restructuring or a
    // different (slower) intersection method, neither of which is really ideal
    let candidates: Vec<usize> = (0..space.len()).collect();
    bron_kerbosch(Vec::new(), candidates, space, radius, cons)
}

// based on
// C. Bron, J. Kerbosch, Finding All Cliques of an Undirected Graph, 1973
// https://doi.org/10.1145/362342.362367
// NOTE O(3^(n / 3)) when max_k = max_dim + 1 = infinity (not accounting for the constructing simplex step)
// NOTE similar to the algorithm described in https://ieeexplore.ieee.org/document/1559964
// TODO try https://arxiv.org/abs/2311.13798v2
// NOTE another idea could be a maximal clique finding algorithm and after breaking down those maximal
// cliques --- Tomita et al. https://doi.org/10.1016/j.tcs.2006.06.015 should be easy enough to implement
// as it uses bron kerbosch as its base
// instead of pivoting at max_(v in P union X) (N(v) intersection (P union X)) as the article
// suggests is optimal a more practical approach might be to just find max_(v in P union X) |N(v)|
// NOTE assuming candidates and adjacency[i] are already sorted
fn bron_kerbosch<M>(
    clique: Vec<usize>,
    candidates: Vec<usize>,
    space: &PointCloud<M>,
    radius: fn(&[usize], &PointCloud<M>, &Construction) -> Option<f64>,
    cons: &mut Construction,
) where
    M: Metric,
{
    if clique.len() > 2 {
        match radius(&clique, space, cons) {
            Some(d) => cons.push(clique.clone(), d),
            None => return,
        }
    }

    if clique.len() == cons.max_dim + 1 || candidates.is_empty() {
        return;
    }
    if clique.len() + candidates.len() < 3 {
        return;
    }

    for (i, &x) in candidates.iter().enumerate() {
        bron_kerbosch(
            join_back(&clique, x),
            intersection_ordered(&candidates[i..], &cons.adjacency[&x]),
            space,
            radius,
            cons,
        )
    }
}

fn intersection_ordered(a: &[usize], b: &[usize]) -> Vec<usize> {
    let (m, n) = (a.len(), b.len());
    let (mut i, mut j) = (0, 0);
    let mut intersection: Vec<usize> = Vec::new();
    while i < m && j < n {
        if a[i] < b[j] {
            i += 1;
            continue;
        }
        if b[j] < a[i] {
            j += 1;
            continue;
        }
        intersection.push(a[i]);
        i += 1;
        j += 1;
    }
    intersection
}

fn join_back(v: &[usize], x: usize) -> Vec<usize> {
    [v, vec![x].as_slice()].concat()
}
