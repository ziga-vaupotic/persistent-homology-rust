use crate::algebra::{ BoundaryMatrices, BoundaryMatrix, ReducedBoundaryMatrix, ReducedBoundaryMatrices };
use std::collections::HashSet;

pub struct PersistencePair {
    pub dimension: usize,
    pub birth: usize,
    pub death: Option<usize>,
}

pub struct PersistenceDiagram {
    pub pairs: Vec<PersistencePair>,
}

pub fn compute_persistence(
    matrix: &ReducedBoundaryMatrix,
    dimension: usize,
) -> Vec<PersistencePair> {
    let mut pairs = Vec::new();

    // Birth simplices that get killed
    let mut paired_births = HashSet::new();

    // Finite intervals
    for (local_col, pivot_row) in matrix.low.iter().enumerate() {
        if let Some(birth) = pivot_row {
            let death = matrix.matrix.column_indices[local_col];

            paired_births.insert(*birth);

            pairs.push(PersistencePair {
                dimension,
                birth: *birth,
                death: Some(death),
            });
        }
    }

    // Infinite intervals:
    // zero reduced columns whose birth simplex
    // never appears as a pivot row.
    for (local_col, col) in matrix.matrix.columns.iter().enumerate() {
        if col.is_empty() {
            let birth = matrix.matrix.column_indices[local_col];

            if !paired_births.contains(&birth) {
                pairs.push(PersistencePair {
                    dimension,
                    birth,
                    death: None,
                });
            }
        }
    }

    pairs
}

pub fn compute_persistence_diagram(
    matrices: &ReducedBoundaryMatrices
) -> PersistenceDiagram {
    let mut pairs = Vec::new();

    for (dimension, matrix) in matrices.iter().enumerate() {
        pairs.extend(compute_persistence(
            matrix,
            dimension,
        ));
    }

    PersistenceDiagram { pairs }
}
