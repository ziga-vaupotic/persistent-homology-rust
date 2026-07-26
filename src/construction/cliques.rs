//! Clique enumeration algorithms for complex construction.
//!
//! This module provides clique finding utilities using the Bron-Kerbosch algorithm
//! to enumerate all maximal cliques in a graph defined by point distances.

use crate::construction::Construction;
use crate::geometry::{Metric, PointCloud};

/// Find all cliques in the adjacency graph defined by the construction state.
///
/// This function uses the Bron-Kerbosch algorithm to enumerate all cliques and applies
/// a user-provided radius function to determine which cliques should be added to the
/// construction. The radius function receives a clique and should return either the
/// filtration value for that clique or `None` if it should be rejected.
///
/// # Arguments
///
/// * `space` - The point cloud with associated metric.
/// * `radius` - Callback function that evaluates whether a clique is valid and returns its filtration value.
/// * `cons` - Mutable reference to the construction state that accumulates valid simplices.
///
/// # Complexity
///
/// $O(3^{n/3})$ in the worst case when enumerating all maximal cliques.
pub fn find_all<const D: usize, M>(
    space: &PointCloud<D, M>,
    radius: fn(&[usize], &PointCloud<D, M>, &Construction) -> Option<f64>,
    cons: &mut Construction,
) where
    M: Metric<D>,
{
    // TODO add degeneracy ordering
    // do not forget that candidates still has to be ordered
    // https://arxiv.org/abs/1006.5440
    // algorithm described in the article runs O(n d 3^(d / 3)) where d is degeneracy
    // https://en.wikipedia.org/wiki/Degeneracy_(graph_theory)#Algorithms
    // NOTE somewhat unviable with the current setup --- would need major restructuring or a
    // different (slower) intersection method, neither of which is really ideal
    bron_kerbosch(
        &mut Vec::new(),
        (0..space.len()).collect(),
        space,
        radius,
        cons,
    )
}

/// Recursive Bron-Kerbosch algorithm for clique enumeration.
///
/// Implementation based on:
/// C. Bron, J. Kerbosch, Finding All Cliques of an Undirected Graph, 1973
/// https://doi.org/10.1145/362342.362367
///
/// # Arguments
///
/// * `clique` - Current clique being extended (initially empty).
/// * `candidates` - Vertices that could extend the current clique.
/// * `space` - The point cloud with associated metric.
/// * `radius` - Callback to evaluate clique validity and filtration value.
/// * `cons` - Mutable construction state.
///
/// # Assumptions
///
/// Assumes both `candidates` and adjacency lists are kept sorted for efficient intersection.
///
/// # Complexity
///
/// $O(3^{n/3})$ when max_k = infinity. Similar to the algorithm described in https://ieeexplore.ieee.org/document/1559964.
///
/// # TODO
///
/// - Implement degeneracy ordering as per https://arxiv.org/abs/1006.5440
/// - Explore pivoting strategies from https://arxiv.org/abs/2311.13798v2
/// - Consider maximal clique finding (Tomita et al.) as a preprocessing step
fn bron_kerbosch<const D: usize, M>(
    clique: &mut Vec<usize>,
    candidates: Vec<usize>,
    space: &PointCloud<D, M>,
    radius: fn(&[usize], &PointCloud<D, M>, &Construction) -> Option<f64>,
    cons: &mut Construction,
) where
    M: Metric<D>,
{
    if clique.len() > 2 {
        match radius(clique, space, cons) {
            Some(d) => cons.push(clique.clone(), d),
            None => return,
        }
    }

    if clique.len() == cons.max_k || candidates.is_empty() || clique.len() + candidates.len() < 3 {
        return;
    }

    for (i, &x) in candidates.iter().enumerate() {
        clique.push(x);
        bron_kerbosch(
            clique,
            intersection_ordered(&candidates[i..], &cons.adjacency[&x]),
            space,
            radius,
            cons,
        );
        clique.pop();
    }
}

/// Compute the intersection of two sorted sequences.
///
/// Returns the set of elements that appear in both sequences.
/// Both input sequences must be sorted in ascending order.
///
/// # Arguments
///
/// * `a` - First sorted sequence.
/// * `b` - Second sorted sequence.
///
/// # Returns
///
/// A new sorted vector containing elements present in both `a` and `b`.
///
/// # Complexity
///
/// $O(m + n)$ where $m$ and $n$ are the lengths of the input sequences.
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
