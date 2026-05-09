
// Gaussian elimination over Z2

pub fn reduce_matrix(matrix: &mut Vec<Vec<usize>>) -> usize {

    let mut max_row = 0;
    for col in matrix.iter() {
        for &r in col {
            if r > max_row {
                max_row = r;
            }
        }
    }
    let num_rows = max_row + 1;
    let mut low: Vec<Option<usize>> = vec![None; num_rows];

    for j in 0..matrix.len() {
        let mut col = matrix[j].clone();
        while !col.is_empty() {
            let pivot_row = col[0];
            if let Some(existing_j) = low[pivot_row] {

                let existing_col = &matrix[existing_j].clone();
                let mut new_col = Vec::new();
                let mut i = 0;
                let mut p = 0;
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
                while i < col.len() {
                    new_col.push(col[i]);
                    i += 1;
                }
                while p < existing_col.len() {
                    new_col.push(existing_col[p]);
                    p += 1;
                }
                col = new_col;
            } else {
                low[pivot_row] = Some(j);
                matrix[j] = col;
                break;
            }
        }
    }
    low.iter().filter(|&&l| l.is_some()).count()

}