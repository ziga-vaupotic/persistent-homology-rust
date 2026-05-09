use crate::algebra::matrices::{BoundaryMatrix, BoundaryMatrices};

pub struct PersistencePair {
    pub birth: usize,
    pub death: Option<usize>,
}

pub struct PersistenceDiagram {
    pub pairs: Vec<PersistencePair>,
}


pub fn compute_persistence_diagram(
    matrix: &BoundaryMatrix,
    low: &[Option<usize>],
) -> PersistenceDiagram {
    let num_cols = matrix.columns().len();

    let mut pairs = Vec::new();
    let mut is_birth = vec![false; num_cols];

    // births → deaths
    for (row, &col) in low.iter().enumerate() {
        if let Some(birth_col) = col {
            is_birth[birth_col] = true;

            pairs.push(PersistencePair {
                birth: birth_col,
                death: Some(row),
            });
        }
    }

    // infinite classes
    for col in 0..num_cols {
        if !is_birth[col] {
            pairs.push(PersistencePair {
                birth: col,
                death: None,
            });
        }
    }

    PersistenceDiagram { pairs }
}