use std::collections::HashMap;

// The boundary matrix is an element of GF(2)
// Currently the columns act as a a way to stor 1s
// E.g. [1,2, 3] means that in column id. 1 there are
// 1s at position 1 2 and 3.

pub struct BoundaryMatrix {
    pub columns: Vec<Vec<usize>>,
    pub column_indices: Vec<usize>, // global filtration indices
}
pub struct ReducedBoundaryMatrix {
    pub matrix: BoundaryMatrix,
    pub low: Vec<Option<usize>>,
}

impl BoundaryMatrix {
    pub fn new(mut columns: Vec<Vec<usize>>, column_indices: Vec<usize>) -> Self {
        assert_eq!(columns.len(), column_indices.len());

        for col in &mut columns {
            col.sort_unstable();
            col.dedup();
        }

        Self {
            columns,
            column_indices,
        }
    }

    pub fn columns(&self) -> &Vec<Vec<usize>> {
        &self.columns
    }

    pub fn reduce(&self) -> ReducedBoundaryMatrix {
        let mut matrix = self.columns.clone();

        // low[j] = largest row index in reduced column j
        let mut low: Vec<Option<usize>> = vec![None; matrix.len()];

        // pivot row -> column owning that pivot
        let mut pivot_owner: HashMap<usize, usize> = HashMap::new();

        for j in 0..matrix.len() {
            let mut col = matrix[j].clone();

            loop {
                if col.is_empty() {
                    break;
                }

                // Standard persistence convention:
                // low(j) = largest row index
                let pivot_row = *col.last().unwrap();

                if let Some(&owner) = pivot_owner.get(&pivot_row) {
                    let owner_col = &matrix[owner];

                    // Z2 symmetric difference
                    let mut new_col = Vec::new();
                    let (mut a, mut b) = (0, 0);

                    while a < col.len() && b < owner_col.len() {
                        match col[a].cmp(&owner_col[b]) {
                            std::cmp::Ordering::Less => {
                                new_col.push(col[a]);
                                a += 1;
                            }
                            std::cmp::Ordering::Greater => {
                                new_col.push(owner_col[b]);
                                b += 1;
                            }
                            std::cmp::Ordering::Equal => {
                                a += 1;
                                b += 1;
                            }
                        }
                    }

                    new_col.extend_from_slice(&col[a..]);
                    new_col.extend_from_slice(&owner_col[b..]);

                    col = new_col;
                } else {
                    low[j] = Some(pivot_row);
                    pivot_owner.insert(pivot_row, j);
                    break;
                }
            }

            matrix[j] = col;
        }

        let mut nonzero = 0;
        let mut zero = 0;

        for col in &matrix {
            if col.is_empty() {
                zero += 1;
            } else {
                nonzero += 1;
            }
        }

        println!("zero={} nonzero={}", zero, nonzero);

        ReducedBoundaryMatrix {
            matrix: BoundaryMatrix {
                columns: matrix,
                column_indices: self.column_indices.clone(),
            },
            low,
        }
    }

    // Computes ranks of the largest chain complex
    pub fn rank(&self) -> usize {
        self.reduce()
            .low
            .into_iter()
            .filter(|x| x.is_some())
            .count()
    }
}
pub type BoundaryMatrices = Vec<BoundaryMatrix>;
pub type ReducedBoundaryMatrices = Vec<ReducedBoundaryMatrix>;

#[cfg(test)]
mod tests {
    use super::*;

    fn test_boundary_matrix(columns: Vec<Vec<usize>>) -> BoundaryMatrix {
        let column_indices = (0..columns.len()).collect();
        BoundaryMatrix::new(columns, column_indices)
    }

    #[test]
    fn test_new_sorts_and_dedups() {
        let bm = test_boundary_matrix(vec![vec![3, 1, 1, 2], vec![2, 2, 0]]);

        assert_eq!(bm.columns(), &vec![vec![1, 2, 3], vec![0, 2],]);
    }

    #[test]
    fn test_rank_independent_columns() {
        let bm = test_boundary_matrix(vec![vec![0], vec![1], vec![2]]);

        assert_eq!(bm.rank(), 3);
    }

    #[test]
    fn test_rank_single_dependency() {
        // c3 = c1 + c2 over GF(2)
        let bm = test_boundary_matrix(vec![vec![0, 1], vec![1], vec![0]]);

        assert_eq!(bm.rank(), 2);
    }

    #[test]
    fn test_rank_full_dependency() {
        let bm = test_boundary_matrix(vec![vec![0], vec![0], vec![0]]);

        assert_eq!(bm.rank(), 1);
    }

    #[test]
    fn test_reduce_pivots_consistent_with_rank() {
        let bm = test_boundary_matrix(vec![vec![0], vec![1], vec![0, 1]]);

        let rank_from_low = bm.reduce().low.iter().filter(|x| x.is_some()).count();
        let rank_direct = bm.rank();

        assert_eq!(rank_from_low, rank_direct);
    }

    #[test]
    fn test_pivots_are_valid_indices() {
        let bm = test_boundary_matrix(vec![vec![0], vec![1], vec![2]]);

        let reduce = bm.reduce();

        for pivot in reduce.low.iter().flatten() {
            assert!(*pivot < bm.columns().len());
        }
    }

    #[test]
    fn test_rank_never_exceeds_columns() {
        let bm = test_boundary_matrix(vec![vec![0, 1], vec![1]]);

        assert!(bm.rank() <= bm.columns().len());
    }
}
