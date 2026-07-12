


use rand::seq::SliceRandom;

use crate::geometry::{ Point, PointSet, Ball };
use crate::construction::common::*;


// NOTE in current state can be nummerically unstable for dimensions > 2
// TODO nummerical stability improvements
// E. Welzl, Smallest enclosing disks (balls and ellipsoids), 1991
// https://doi.org/10.1007/BFb0038202
// NOTE should be only used for low dimensions (maybe up to 10?)
// calculates ball from boundary from scratch every time so Gärtner 1999 should be prefered in most
// cases (once its implemented)
pub fn welzl(points : &Vec<usize>, space : &PointSet) -> Ball {
    let mut P = points.clone();
    P.shuffle(&mut rand::rng());

    welzl_rec(&mut P, Vec::new(), space)
}


// using Gärtner 1999 move-to-front implementation
// NOTE should probably be using LinkedList instead of Vec
// unfortunatelly rust LinkedLists do not support remove/insert operations (nightly versions only)
// should be OK as long as |P| does not exceed 10
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


fn from_boundary(boundary : &Vec<usize>, space : &PointSet) -> Ball {
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
// find affine subspace containing boundary find isometry to R^(n - 1), where n is size of boundary
// calculate miniball there then move the center back to original subspace
// no need to change the radius as we have an isometry
pub fn on_affine_subspace(boundary : &Vec<usize>, space : &PointSet) -> Ball {
    let dim = space.dim();

    let linear_parts : Vec<Point> = (1..boundary.len())
        .map(|x| Point::difference(space.get(boundary[x]), space.get(boundary[0]))).collect();

    // find Q such that Q(q_i - q_0) = (a_1, ..., a_n, 0, ..., 0)
    let (Q, LI) = extend_to_basis(&linear_parts); // orthonormal matrix => isometry
    let n = LI.len();

    let mut new_space_points : Vec<Point> = Vec::new();
    for i in LI {
        let mut u = multiply(&Q, &linear_parts[i].coords);
        u.truncate(n); // all components after are zero as per Gram Schmidt
        new_space_points.push(Point::new(u));
    }
    new_space_points.push(Point::new(vec![0.0; n])); // add the n + 1 th point
    let new_space = PointSet::new_no_check(new_space_points);
    let new_boundary : Vec<usize> = (0..(n + 1)).collect();

    // find center in subpace spanned by {Q(q_i - q_0)}_i
    let miniball = circumsphere(&new_boundary, &new_space);

    // get back the original center with Q^T center_new + q_0 = center
    let mut center_new = miniball.o().coords.clone();
    let mut zeros = vec![0.0; dim - n];
    center_new.append(&mut zeros);
    let mut center = Point::new(multiply_transpose(&Q, &center_new));
    center.add(space.get(boundary[0]));

    Ball::new(center, miniball.radius)
}


fn multiply(A : &Vec<Vec<f64>>, v : &Vec<f64>) -> Vec<f64> {
    let m = A.len();
    let n = A[0].len();

    let mut result : Vec<f64> = Vec::new();
    for i in 0..m {
        let s = (0..n).map(|j| v[j] * A[i][j]).sum::<f64>();
        result.push(s);
    }
    result
}


fn multiply_transpose(A : &Vec<Vec<f64>>, v : &Vec<f64>) -> Vec<f64> {
    let m = A.len();
    let n = A[0].len();

    let mut result : Vec<f64> = Vec::new();
    for i in 0..n {
        let s = (0..m).map(|j| v[j] * A[j][i]).sum::<f64>();
        result.push(s);
    }
    result
}


// https://en.wikipedia.org/wiki/Gram%E2%80%93Schmidt_process
// NOTE using modified Gram Schmidt process --- might still be nummerically unstable
pub fn extend_to_basis(points : &Vec<Point>) -> (Vec<Vec<f64>>, Vec<usize>) {
    let dim = points[0].dim();
    let n = points.len();

    let mut base : Vec<Point> = Vec::new();

    let mut u_0 = points[0].clone();
    u_0.normalize();
    base.push(u_0);

    for i in 1..n {
        base.push(points[i].clone());
    }
    for i in 0..dim { // extend to span of R^dim
        base.push(Point::standard_unit(i, dim));
    }

    // Gram Schmidt
    for i in 1..base.len() {
        for j in i..base.len() {
            let projection = Point::projection_normal(&base[j], &base[i - 1]);
            base[j].subtract(&projection);
        }
        base[i].normalize();
    }

    let linearly_independent : Vec<usize> = (0..n).filter(|&x| !base[x].is_zero()).collect();
    (base.into_iter().filter(|x| !x.is_zero()).map(|x| x.coords.clone()).collect(), linearly_independent)
}


// generalised formula from https://mathworld.wolfram.com/Circumsphere.html
// as well as https://en.wikipedia.org/wiki/Circumcircle both taken on the 11th of July 2026
fn circumsphere(boundary : &Vec<usize>, space : &PointSet) -> Ball {
    let dim = space.dim();
    let n = dim + 1; // length of boundary, if longer just take the first n

    fn det(A : &mut Vec<Vec<f64>>) -> f64 {
        if A.len() < 6 { return det_naive(A); }
        det_LU(A)
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


// n! < n^3 for n < 6 => faster than LU for small n
fn det_naive(A : &mut Vec<Vec<f64>>) -> f64 { // mut just to match with det_LU
    if A.len() == 1 { return A[0][0]; }
    if A.len() == 2 { return A[0][0] * A[1][1] - A[0][1] * A[1][0]; }
    if A.len() == 3 {
        return A[0][0] * (A[1][1] * A[2][2] - A[2][1] * A[1][2])
            - A[0][1] * (A[1][0] * A[2][2] - A[2][0] * A[1][2])
            + A[0][2] * (A[1][0] * A[2][1] - A[2][0] * A[1][1]);
    }

    let mut result = 0.0;
    for i in 0..A.len() {
        let mut A_i = A.clone();
        A_i.remove(0);
        for j in 0..A_i.len() {
            A_i[j].remove(i);
        }
        result += (-1.0_f64).powf(i as f64) * A[0][i] * det_naive(&mut A_i);
    }
    result
}


fn det_LU(mut A : &mut Vec<Vec<f64>>) -> f64 {
    let n = A.len();
    let mut U = vec![vec![0.0; n]; n]; // only interested in the values on the diagonal
    let mut L = vec![vec![0.0; n]; n];

    LU_decomposition(&mut L, &mut U, &mut A);

    let mut det = 1.0;
    for i in 0..n {
        det *= L[i][i];
    }
    det
}


// assumes that the appropriate values of U and L ie. upper and lower triangles
// respectively are already zero before input (not a problem if not)
fn LU_decomposition(L : &mut Vec<Vec<f64>>, U : &mut Vec<Vec<f64>>, A : &mut Vec<Vec<f64>>) {
    let n = A.len();
    for i in 0..n {
        for j in i..n {
            L[j][i] = A[j][i];
            for k in 0..i {
                L[j][i] = L[j][i] - L[j][k] * U[k][i];
            }
        }
        for j in i..n {
            if L[i][i].abs() < 1e-12 { // tolerance
                U[i][j] = 0.0;
                continue;
            }
            U[i][j] = A[i][j] / L[i][i];
            for k in 0..i {
                U[i][j] = U[i][j] - ((L[i][k] * U[k][j]) / L[i][i]);
            }
        }
    }
}
