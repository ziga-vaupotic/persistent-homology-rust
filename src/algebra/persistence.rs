use crate::algebra::matrices::{BoundaryMatrices, BoundaryMatrix, ReducedBoundaryMatrix, ReducedBoundaryMatrices};
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
    let num_cols = matrix.matrix.columns().len();

    let mut pairs = Vec::new();
    let mut paired_births = vec![false; num_cols];

    // Finite intervals
    for (death_col, &birth_col) in matrix.low.iter().enumerate() {
        if let Some(birth_col) = birth_col {
            paired_births[birth_col] = true;

            pairs.push(PersistencePair {
                dimension,
                birth: birth_col,
                death: Some(death_col),
            });
        }
    }

    // Infinite intervals
    for birth_col in 0..num_cols {
        if !paired_births[birth_col] {
            pairs.push(PersistencePair {
                dimension,
                birth: birth_col,
                death: None,
            });
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