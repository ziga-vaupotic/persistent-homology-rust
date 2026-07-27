//! Discretises filtration simplices into algebraic boundary matrices.
//!
//! Each simplex boundary is converted into a sparse column representation of a
//! boundary operator over `\mathbb{Z}_2`.

use crate::topology::{Cell, Filtration, Simplex};

use crate::algebra::matrices::{BoundaryMatrices, BoundaryMatrix, ReducedBoundaryMatrices};

use std::collections::HashMap;

/// Build the boundary matrices for a filtration.
///
/// Each returned `BoundaryMatrix` corresponds to a dimension `d` and contains columns
/// representing the boundary of each `d`-simplex. The `column_indices` preserve the
/// original global filtration order of simplices.
/// Build the boundary matrices for a filtration.
///
/// Creates a sparse boundary matrix for each dimension $d \geq 1$. Each column represents
/// the boundary of a $d$-simplex over $\mathbb{Z}_2$ (i.e., modulo 2), stored as a set of
/// row indices where the matrix has 1 entries. The `column_indices` in each matrix preserve
/// the original global filtration ordering of simplices for later persistence computation.
///
/// # Arguments
///
/// * `filtration` - The ordered filtration of simplices.
///
/// # Returns
///
/// A vector of `BoundaryMatrix` structures, indexed by dimension. Empty dimensions are skipped.
///
/// # Example
///
/// ```ignore
/// use persistent_homology::algebra::discretisation::build_boundary_matrices;
/// use persistent_homology::topology::{Simplex, Filtration};
///
/// let simplices = vec![
///     Simplex::new(vec![0], 0.0),
///     Simplex::new(vec![1], 0.0),
///     Simplex::new(vec![0, 1], 1.0),
/// ];
/// let filtration = Filtration::new(simplices);
/// let matrices = build_boundary_matrices(&filtration);
/// ```
pub fn build_boundary_matrices(filtration: &Filtration<Simplex>) -> BoundaryMatrices {
    let mut by_dim: Vec<Vec<&Simplex>> = Vec::new();

    for simplex in &filtration.cells {
        let d = simplex.dim();

        if d >= by_dim.len() {
            by_dim.resize(d + 1, Vec::new());
        }

        by_dim[d].push(simplex);
    }

    // Global filtration index of every simplex
    let simplex_to_index: HashMap<&Simplex, usize> = filtration
        .cells
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
/// Reduce each boundary matrix to compute pivots for persistence computation.
///
/// Applies the standard persistent homology reduction algorithm to each matrix,
/// computing the (low) pivot indices required for extracting persistence pairs.
/// The reduction occurs over $\mathbb{Z}_2$ using Gaussian elimination with column operations.
///
/// # Arguments
///
/// * `matrices` - A vector of unreduced boundary matrices.
///
/// # Returns
///
/// A vector of reduced boundary matrices, one for each input matrix, with pivot indices computed.
///
/// # Example
///
/// ```ignore
/// use persistent_homology::algebra::discretisation::{build_boundary_matrices, reduce_boundary_matrices};
/// use persistent_homology::topology::{Simplex, Filtration};
///
/// let filtration = Filtration::new(vec![
///     Simplex::new(vec![0], 0.0),
///     Simplex::new(vec![1], 0.5),
///     Simplex::new(vec![0, 1], 1.0),
/// ]);
/// let matrices = build_boundary_matrices(&filtration);
/// let reduced = reduce_boundary_matrices(&matrices);
/// ```
pub fn reduce_boundary_matrices(matrices: &BoundaryMatrices) -> ReducedBoundaryMatrices {
    matrices.iter().map(|matrix| matrix.reduce()).collect()
}
