use crate::topology::filtration::Filtration;
use crate::topology::simplex::Simplex;

use crate::algebra::types::{BoundaryMatrix, BoundaryMatrices};


use std::collections::HashMap;

// A chain complex induces a family of boundary operators \partial_j: C_j \to C_{j-1}
// We discretise them each one of them with a boundary matrix.


pub fn build_boundary_matrices(
    filtration: &Filtration,
) -> BoundaryMatrices {
    let mut by_dim: Vec<Vec<&Simplex>> = Vec::new();

    for simplex in &filtration.simplices {
        let d = simplex.dimension();
        if d >= by_dim.len() {
            by_dim.resize(d + 1, Vec::new());
        }
        by_dim[d].push(simplex);
    }


    // to each simplex we attach a hash to index them

    let mut simplex_to_index: HashMap<&Simplex, usize> = HashMap::new();
    for (i, simplex) in filtration.simplices.iter().enumerate() {
        simplex_to_index.insert(simplex, i);
    }

    let mut matrices: BoundaryMatrices = Vec::new();

    for k in 1..by_dim.len() {

        // boundary matrix takes d-dimensioanl simplexes to d-1 dim simplexes.

        let cols = &by_dim[k];

        let mut columns: Vec<Vec<usize>> = vec![Vec::new(); cols.len()];

        for (j, &simplex) in cols.iter().enumerate() {
            for (_, sub) in simplex.boundary() {
                if let Some(&row_idx) = simplex_to_index.get(&sub) {
                    columns[j].push(row_idx);
                }
            }
            columns[j].sort();
        }

        matrices.push(BoundaryMatrix::new(columns));
    }

    matrices
}