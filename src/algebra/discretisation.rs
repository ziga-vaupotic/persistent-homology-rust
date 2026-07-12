use crate::topology::{ Filtration, Simplex };

use crate::algebra::{ BoundaryMatrices, BoundaryMatrix, ReducedBoundaryMatrices, ReducedBoundaryMatrix };

use std::collections::HashMap;

// A chain complex induces a family of boundary operators \partial_j: C_j \to C_{j-1}
// We discretise them each one of them with a boundary matrix.


pub fn build_boundary_matrices(
    filtration: &Filtration,
) -> BoundaryMatrices {
    let mut by_dim: Vec<Vec<&Simplex>> = Vec::new();

    for simplex in &filtration.simplices {
        let d = simplex.dim();

        if d >= by_dim.len() {
            by_dim.resize(d + 1, Vec::new());
        }

        by_dim[d].push(simplex);
    }

    // Global filtration index of every simplex
    let simplex_to_index: HashMap<&Simplex, usize> =
        filtration
            .simplices
            .iter()
            .enumerate()
            .map(|(i, s)| (s, i))
            .collect();

    let mut matrices = Vec::new();

    for k in 1..by_dim.len() {
        let cols = &by_dim[k];

        let mut columns = Vec::with_capacity(cols.len());
        let mut column_indices = Vec::with_capacity(cols.len());

        for simplex in cols {
            let mut column = Vec::new();

            for (_, face) in simplex.boundary() {
                if let Some(&row_idx) = simplex_to_index.get(&face) {
                    column.push(row_idx);
                }
            }

            column.sort();

            columns.push(column);

            // Store the GLOBAL filtration index
            column_indices.push(simplex_to_index[simplex]);
        }

        matrices.push(
            BoundaryMatrix {
                columns,
                column_indices,
            }
        );
    }

    matrices
}


pub fn reduce_boundary_matrices(
    matrices: &BoundaryMatrices,
) -> ReducedBoundaryMatrices {
    matrices
        .iter()
        .map(|matrix| {
            matrix.reduce()
        })
        .collect()
}

