//! Discretises filtration simplices into algebraic boundary matrices.
//!
//! Each simplex boundary is converted into a sparse column representation of a
//! boundary operator over `\mathbb{Z}_2`.

use crate::topology::{Filtration, Simplex};

use crate::algebra::matrices::{BoundaryMatrices, BoundaryMatrix, ReducedBoundaryMatrices};

use std::collections::HashMap;

/// Build the boundary matrices for a filtration.
///
/// Each returned `BoundaryMatrix` corresponds to a dimension `d` and contains columns
/// representing the boundary of each `d`-simplex. The `column_indices` preserve the
/// original global filtration order of simplices.
pub fn build_boundary_matrices(filtration: &Filtration) -> BoundaryMatrices {
    let mut by_dim: Vec<Vec<&Simplex>> = Vec::new();

    for simplex in &filtration.simplices {
        let d = simplex.dim();

        if d >= by_dim.len() {
            by_dim.resize(d + 1, Vec::new());
        }

        by_dim[d].push(simplex);
    }

    // Global filtration index of every simplex
    let simplex_to_index: HashMap<&Simplex, usize> = filtration
        .simplices
        .iter()
        .enumerate()
        .map(|(i, s)| (s, i))
        .collect();

    let mut matrices = Vec::new();

    for cols in by_dim.iter().skip(1) {
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

        matrices.push(BoundaryMatrix {
            columns,
            column_indices,
        });
    }

    matrices
}

/// Reduce each boundary matrix to its canonical `\mathbb{Z}_2` persistence form.
///
/// The reduction uses the standard persistence convention, where the pivot (low) entry
/// of each nonzero column is the largest row index present in that column.
pub fn reduce_boundary_matrices(matrices: &BoundaryMatrices) -> ReducedBoundaryMatrices {
    matrices.iter().map(|matrix| matrix.reduce()).collect()
}
