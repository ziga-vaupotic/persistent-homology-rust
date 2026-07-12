


// bad practice to have a utils file? --- might change at some point


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


pub fn join_back(v : &Vec<usize>, x : usize) -> Vec<usize> {
    [v, vec![x].as_slice()].concat()
}


pub fn move_front(v : &mut Vec<usize>, x : usize) {
    let u = v[x];
    v.remove(x);
    v.splice(..0, [u]);
}


pub fn cut_at(P : &Vec<usize>, x : usize) -> Vec<usize> {
    let mut P_new = P.clone();
    P_new.truncate(x);
    P_new
}
