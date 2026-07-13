


// bad practice having a utils file? might change at some point


pub fn intersection_ordered(a : &Vec<usize>, b : &Vec<usize>) -> Vec<usize> {
    let (m, n) = (a.len(), b.len());
    let (mut i, mut j) = (0, 0);
    let mut intersection = vec![];
    while i < m && j < n {
        if a[i] < b[j] {
            i += 1;
            continue
        }
        if b[j] < a[i] {
            j += 1;
            continue
        }
        intersection.push(a[i]);
        i += 1;
        j += 1;
    }
    intersection
}
