use crate::geometry::point::Point;
use crate::geometry::point_set::PointSet;

use crate::topology::simplex::Simplex;
use crate::topology::filtration::Filtration;

use std::collections::{HashMap, HashSet};

pub fn vietoris_rips_flag(points: &PointSet, max_dim: usize) -> Filtration {
    let n = points.len();
    let point_dim = points.dim();

    let mut edges: Vec<(usize, usize, f64)> = Vec::new();

    for i in 0..n {
        for j in (i + 1)..n {
            let d = if point_dim > 0 {
                points.get(i).distance(&points.get(j))
            } else {
                points.get(i).distance(&points.get(j))
            };
            edges.push((i, j, d));
        }
    }

    edges.sort_by(|a, b| a.2.partial_cmp(&b.2).unwrap());

    // Build adjacency (grows with filtration)
    let mut adj: Vec<HashSet<usize>> = vec![HashSet::new(); n];

    let mut simplices: Vec<Simplex> = Vec::new();

    // All vertices (dimension 0)
    for i in 0..n {
        simplices.push(Simplex::new(vec![i], 0.0));
    }

    // Map clique to filtration value
    let mut simplex_map: HashMap<Vec<usize>, f64> = HashMap::new();

    for (i, j, w) in edges {
        adj[i].insert(j);
        adj[j].insert(i);

        let edge = vec![i.min(j), i.max(j)];
        simplices.push(Simplex::new(edge.clone(), w));
        simplex_map.insert(edge, w);

        // Try to extend cliques up to max_dim
        let mut new_cliques = vec![vec![i, j]];

        while let Some(clique) = new_cliques.pop() {
            if clique.len() > max_dim + 1 {
                continue;
            }

            let mut common_neighbors: HashSet<usize> = adj[clique[0]].clone();

            for &v in &clique[1..] {
                common_neighbors = common_neighbors
                    .intersection(&adj[v])
                    .cloned()
                    .collect();
            }

            for &v in &common_neighbors {
                if clique.contains(&v) {
                    continue;
                }

                let mut new_clique = clique.clone();
                new_clique.push(v);
                new_clique.sort_unstable();
                new_clique.dedup();

                if new_clique.len() <= max_dim + 1 {
                    let filtration = w;

                    simplices.push(Simplex::new(new_clique.clone(), filtration));
                    new_cliques.push(new_clique);
                }
            }
        }
    }


    simplices.sort_by(|a, b| {
        a.filtration_value
            .partial_cmp(&b.filtration_value)
            .unwrap()
            .then(a.dimension().cmp(&b.dimension()))
            .then(a.vertices.cmp(&b.vertices))
    });

    Filtration { simplices }
}