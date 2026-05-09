
// The boundary matrix is an element of GF(2)
// Currently the columns act as a a way to stor 1s
// E.g. [1,2, 3] means that in column id. 1 there are
// 1s at position 1 2 and 3.

pub struct BoundaryMatrix {
    columns: Vec<Vec<usize>>,
}

impl BoundaryMatrix {
    pub fn new(mut columns: Vec<Vec<usize>>) -> Self {
        // enforce consistency i.e. row-ids are sorted

        for col in &mut columns {
            col.sort_unstable();
            col.dedup();
        }
        Self { columns }
    }

    pub fn columns(&self) -> &Vec<Vec<usize>> {
        &self.columns
    }

    pub fn reduce(&self) -> (BoundaryMatrix, Vec<Option<usize>>) {
        let mut matrix = self.columns.clone();

        let mut max_row = 0;
        for col in &matrix {
            for &r in col {
                max_row = max_row.max(r);
            }
        }

        let num_rows = max_row + 1;
        let mut low: Vec<Option<usize>> = vec![None; num_rows];

        for j in 0..matrix.len() {
            let mut col = matrix[j].clone();

            while !col.is_empty() {
                let pivot_row = col[0];

                if let Some(existing_j) = low[pivot_row] {
                    let existing_col = matrix[existing_j].clone();

                    let mut new_col = Vec::new();
                    let (mut i, mut p) = (0, 0);

                    while i < col.len() && p < existing_col.len() {
                        if col[i] < existing_col[p] {
                            new_col.push(col[i]);
                            i += 1;
                        } else if col[i] > existing_col[p] {
                            new_col.push(existing_col[p]);
                            p += 1;
                        } else {
                            i += 1;
                            p += 1;
                        }
                    }

                    new_col.extend_from_slice(&col[i..]);
                    new_col.extend_from_slice(&existing_col[p..]);

                    col = new_col;
                } else {
                    low[pivot_row] = Some(j);
                    matrix[j] = col;
                    break;
                }
            }
        }

        // return pivots aswell
        (BoundaryMatrix { columns: matrix }, low)
    }


    // Computes ranks of the largest chain complex
    pub fn rank(&self) -> usize {
        let (_, low) = self.reduce();
        low.into_iter().filter(|x| x.is_some()).count()
    }

}
pub type BoundaryMatrices = Vec<BoundaryMatrix>;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_new_sorts_and_dedups() {
        let bm = BoundaryMatrix::new(vec![
            vec![3, 1, 1, 2],
            vec![2, 2, 0],
        ]);

        assert_eq!(bm.columns(), &vec![
            vec![1, 2, 3],
            vec![0, 2],
        ]);
    }

    #[test]
    fn test_rank_independent_columns() {
        let bm = BoundaryMatrix::new(vec![
            vec![0],
            vec![1],
            vec![2],
        ]);

        assert_eq!(bm.rank(), 3);
    }

    #[test]
    fn test_rank_single_dependency() {
        // c3 = c1 + c2 over GF(2)
        let bm = BoundaryMatrix::new(vec![
            vec![0, 1],
            vec![1],
            vec![0],
        ]);

        assert_eq!(bm.rank(), 2);
    }

    #[test]
    fn test_rank_full_dependency() {
        let bm = BoundaryMatrix::new(vec![
            vec![0],
            vec![0],
            vec![0],
        ]);

        assert_eq!(bm.rank(), 1);
    }

    #[test]
    fn test_reduce_pivots_consistent_with_rank() {
        let bm = BoundaryMatrix::new(vec![
            vec![0],
            vec![1],
            vec![0, 1],
        ]);

        let (_, low) = bm.reduce();

        let rank_from_low = low.iter().filter(|x| x.is_some()).count();
        let rank_direct = bm.rank();

        assert_eq!(rank_from_low, rank_direct);
    }

    #[test]
    fn test_pivots_are_valid_indices() {
        let bm = BoundaryMatrix::new(vec![
            vec![0],
            vec![1],
            vec![2],
        ]);

        let (_, low) = bm.reduce();

        for pivot in low.iter().flatten() {
            assert!(*pivot < bm.columns().len());
        }
    }

    #[test]
    fn test_rank_never_exceeds_columns() {
        let bm = BoundaryMatrix::new(vec![
            vec![0, 1],
            vec![1],
        ]);

        assert!(bm.rank() <= bm.columns().len());
    }
}