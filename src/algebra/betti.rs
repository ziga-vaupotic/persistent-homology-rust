use crate::algebra::boundary_matrix::build_boundary_matrices;
use crate::algebra::reduction::reduce_matrix;
use crate::topology::Filtration;

pub fn compute_betti_numbers(filtration: &Filtration) -> Vec<usize> {
    let mut matrices = build_boundary_matrices(filtration);

    let mut ranks = vec![0; matrices.len() + 1];

    for (i, matrix) in matrices.iter_mut().enumerate() {
        ranks[i + 1] = reduce_matrix(matrix);
    }

    let max_dim = filtration
        .simplices
        .iter()
        .map(|s| s.dimension())
        .max()
        .unwrap_or(0);

    let mut counts = vec![0; max_dim + 1];
    for s in &filtration.simplices {
        counts[s.dimension()] += 1;
    }

    let mut betti = vec![0; max_dim + 1];

    if max_dim >= 0 {
        betti[0] = counts[0] - ranks[1];
    }

    for k in 1..=max_dim {
        if k + 1 < ranks.len() {
            betti[k] = ranks[k] - ranks[k + 1];
        } else {
            betti[k] = ranks[k];
        }
    }

    betti
}