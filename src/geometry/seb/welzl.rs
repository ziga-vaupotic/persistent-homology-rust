


use nalgebra::{DMatrix, DVector};
use rand::seq::SliceRandom;

use crate::geometry::{ Point, PointSet, Ball };


// E. Welzl, Smallest enclosing disks (balls and ellipsoids), 1991
// https://doi.org/10.1007/BFb0038202
// NOTE in current state can be numerically unstable for dimensions > 2 because of matrix operations
// NOTE should be only used for low dimensions ie. < 5 (maybe up to 10? with a very small space)
// the basic algorith listed in the article is O((dim + 1)(dim + 1)!n) without accounting for the ball
// from boundary step, while time complexity for this specific version of the algorithm is not
// provided, it cannot be much better
// this specific implementation also calculates the ball from a given boundary from scratch every
// time so Gärtner 1999 should be prefered in most cases (once its implemented)
pub fn welzl(points : &Vec<usize>, space : &PointSet) -> Ball {
    let mut P = points.clone();
    P.shuffle(&mut rand::rng());

    welzl_rec(&mut P, Vec::new(), space)
}


// using Gärtner 1999 move-to-front implementation
// NOTE should probably be using LinkedList instead of Vec
// unfortunatelly rust LinkedLists do not support remove/insert operations (nightly versions only)
fn welzl_rec(P : &mut Vec<usize>, B : Vec<usize>, space : &PointSet) -> Ball {
    let mut miniball = from_boundary(&B, space);
    if P.is_empty() || B.len() == space.dim() + 1 {
        return miniball;
    }

    let n = P.len();
    for i in 0..n {
        if miniball.contains(space.get(P[i])) { continue; }
        miniball = welzl_rec(&mut cut_at(P, i), join_back(&B, P[i]), space);
        move_front(P, i);
    }
    miniball
}


fn join_back(v : &Vec<usize>, x : usize) -> Vec<usize> {
    [v, vec![x].as_slice()].concat()
}


fn move_front(v : &mut Vec<usize>, x : usize) {
    let u = v[x];
    v.remove(x);
    v.splice(..0, [u]);
}


fn cut_at(P : &Vec<usize>, x : usize) -> Vec<usize> {
    let mut P_new = P.clone();
    P_new.truncate(x);
    P_new
}


fn from_boundary(boundary : &Vec<usize>, space : &PointSet) -> Ball { // |boundary| <= dim + 1
    match boundary.len() {
        0 => return Ball::new(Point::new(Vec::new()), 0.0),
        1 => return Ball::new(space.get(boundary[0]).clone(), 0.0),
        2 => {
            let mut o = space.sum(boundary);
            o.multiply(1.0 / 2.0);
            let r = o.distance(space.get(boundary[0]));

            return Ball::new(o, r);
       },
       _ => {}
    }

    if boundary.len() < space.dim() + 1 {
        return on_affine_subspace(boundary, space);
    }

    circumsphere(boundary, space)
}


// procedure :
// find smallest affine subspace containing boundary
// find isometry to R^n subset R^d, where n is the dimension of that subspace
// calculate miniball there then move the center back to original subspace
// no need to change the radius as we have an isometry
fn on_affine_subspace(boundary : &Vec<usize>, space : &PointSet) -> Ball {
    let q0 = space.get(boundary[0]);
    let dim = q0.dim();

    let linear_parts = DMatrix::from_fn(
        dim,
        boundary.len() - 1,
        |row, col| {
            space.get(boundary[col + 1]).coords[row]
                - q0.coords[row]
        },
    );

    let (basis, n) = extend_to_basis(&linear_parts);

    // Transform the vectors into the basis coordinates.
    let transformed = basis.transpose() * &linear_parts;

    let mut new_space_points: Vec<Point> = (0..n.min(linear_parts.ncols()))
        .map(|i| {
            let coords: Vec<f64> = (0..n).map(|j| transformed[(j, i)]).collect();
            Point::new(coords)
        })
        .collect();
    new_space_points.push(Point::new(vec![0.0; n]));

    let new_space = PointSet::new_no_check(new_space_points);
    let new_boundary : Vec<usize> = (0..(n + 1)).collect();

    // find center in subpace spanned by {Q(q_i - q_0)}_i
    let miniball = circumsphere(&new_boundary, &new_space);

    // center = basis * center_new + q_0
    let center_new_vec = DVector::from_vec(miniball.o().coords.clone());
    let center_in_original_space = basis * center_new_vec;
    let mut center = Point::new(center_in_original_space.iter().copied().collect());
    center.add(q0);

    Ball::new(center, miniball.radius)
}


fn extend_to_basis(points: &DMatrix<f64>) -> (DMatrix<f64>, usize) {
    let dim = points.nrows();
    let n = points.ncols();

    let matrix = DMatrix::from_fn(dim, n, |row, col| points[(row, col)]);
    let qr = matrix.qr();
    let q = qr.q();
    let rank = numerical_rank(&qr.r(), 1e-12);

    let basis = q.columns(0, rank).into_owned();
    (basis, rank)
}

fn numerical_rank(r: &DMatrix<f64>, tolerance: f64) -> usize {
    let diagonal_len = r.nrows().min(r.ncols());

    let max_diag = (0..diagonal_len)
        .map(|i| r[(i, i)].abs())
        .fold(0.0_f64, f64::max);

    if max_diag == 0.0 {
        return 0;
    }

    (0..diagonal_len)
        .filter(|&i| r[(i, i)].abs() > tolerance * max_diag)
        .count()
}


// generalised formula from https://mathworld.wolfram.com/Circumsphere.html
// https://en.wikipedia.org/wiki/Circumcircle
fn circumsphere(boundary : &Vec<usize>, space : &PointSet) -> Ball {
    let dim = space.dim();
    let n = dim + 1; // length of boundary

    fn det(A : &mut Vec<Vec<f64>>) -> f64 {
        let n = A.len();
        if n == 0 { return 1.0; }
        let matrix = DMatrix::from_fn(n, n, |row, col| A[row][col]);
        matrix.determinant()
    }

    let mut transpose : Vec<Vec<f64>> = Vec::new();
    for i in 0..dim {
        let x_js : Vec<f64> = (0..n).map(|x| space.get(boundary[x]).coords[i]).collect();
        transpose.push(x_js);
    }
    let norms : Vec<f64> = (0..n).map(|x| space.get(boundary[x]).norm_square()).collect();
    let ones = vec![1.0; n];

    let mut c : Vec<f64> = Vec::new();
    for i in 0..dim {
        let mut M_i : Vec<Vec<f64>> = Vec::new(); // D_x, D_y, ... in wolfram reference
        M_i.push(norms.clone());
        for j in 0..dim {
            if j == i { continue; }
            M_i.push(transpose[j].clone());
        }
        M_i.push(ones.clone());
        c.push((-1.0_f64).powf(i as f64) * det(&mut M_i));
    }

    let mut A : Vec<Vec<f64>> = (0..dim).map(|x| transpose[x].clone()).collect();
    A.push(ones);
    let a = 1.0 / (2.0 * det(&mut A));

    let mut center = Point::new(c);
    center.multiply(a);

    let radius = center.distance(space.get(boundary[0]));

    Ball::new(center, radius)
}




#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn circumsphere_for_right_triangle_matches_expected() {
        let points = vec![
            Point::new(vec![0.0, 0.0]),
            Point::new(vec![1.0, 0.0]),
            Point::new(vec![0.0, 1.0]),
        ];
        let space = PointSet::new_no_check(points);
        let boundary = vec![0, 1, 2];

        let ball = circumsphere(&boundary, &space);

        assert!((ball.o().coords[0] - 0.5).abs() < 1e-10);
        assert!((ball.o().coords[1] - 0.5).abs() < 1e-10);
        assert!((ball.radius - (2.0_f64).sqrt() / 2.0).abs() < 1e-10);
    }
}
