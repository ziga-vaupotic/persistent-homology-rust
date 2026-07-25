


use nalgebra::DMatrix;
use rand::seq::SliceRandom;

use crate::geometry::{ Point, PointCloud, Ball, Euclidean };


// E. Welzl, Smallest enclosing disks (balls and ellipsoids), 1991
// https://doi.org/10.1007/BFb0038202
// NOTE in current state can be numerically unstable for dimensions > 2 because of matrix operations
// NOTE should be only used for low dimensions ie. < 5 (maybe up to 10? with a very small space)
// the basic algorith listed in the article is O((dim + 1)(dim + 1)!n) without accounting for the ball
// from boundary step, while time complexity for this specific version of the algorithm is not
// provided, it cannot be much better
// this specific implementation also calculates the ball from a given boundary from scratch every
// time so Gärtner 1999 should be prefered in most cases (once its implemented)
pub fn welzl<M> (points : &Vec<usize>, space : &PointCloud<M>) -> Ball 
where 
    M: Euclidean
{
    let mut P = points.clone();
    P.shuffle(&mut rand::rng());

    welzl_rec(&mut P, Vec::new(), space)
}


// using Gärtner 1999 move-to-front implementation
// NOTE should probably be using LinkedList instead of Vec
// unfortunatelly rust LinkedLists do not support remove/insert operations (nightly versions only)
fn welzl_rec<M> (P : &mut Vec<usize>, B : Vec<usize>, space : &PointCloud<M>) -> Ball
where 
    M: Euclidean
{
    let mut miniball = from_boundary(&B, space);
    let p = miniball.clone();
    println!("{:?} {:?}", p.o().coords, p.r());
    if P.is_empty() || B.len() == space.dim() + 1 {
        return miniball;
    }

    let n = P.len();
    for i in 0..n {
        if space.contained_in_ball(&miniball, space.get(P[i])) { continue; }
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


fn from_boundary<M> (boundary : &Vec<usize>, space : &PointCloud<M>) -> Ball // |boundary| <= dim + 1
where 
    M: Euclidean
{
    match boundary.len() {
        0 => return Ball::new(Point::new(Vec::new()), 0.0),
        1 => return Ball::new(space.get(boundary[0]).clone(), 0.0),
        2 => {
            let mut o = space.get(boundary[0]) + space.get(boundary[1]);
            o.multiply(1.0 / 2.0);
            let r = space.distance(&o, space.get(boundary[0]));

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
fn on_affine_subspace<M> (boundary : &Vec<usize>, space : &PointCloud<M>) -> Ball
where 
    M: Euclidean
{
    let q0 = space.get(boundary[0]);
    let dim = q0.len();

    let linear_parts = DMatrix::from_fn(
        dim,
        boundary.len() - 1,
        |row, col| {
            space.get(boundary[col + 1]).coords[row]
                - q0.coords[row]
        },
    );

    let basis = extend_to_basis(&linear_parts);
    let n = basis.ncols();

    let transformed = basis.transpose() * &linear_parts;

    let mut new_space_points: Vec<Point> = (0..n)
        .map(|i| {
            let coords: Vec<f64> = (0..n).map(|j| transformed[(j, i)]).collect();
            Point::new(coords)
        })
        .collect();
    new_space_points.push(Point::new(vec![0.0; n]));

    let new_space = PointCloud::new_no_check(new_space_points, space.get_geometry());
    let new_boundary : Vec<usize> = (0..(n + 1)).collect();

    // find center in subpace spanned by {Q(q_i - q_0)}_i
    let miniball = circumsphere(&new_boundary, &new_space);

    // center = basis * center_new + q_0
    let center_new_vec = miniball.o().coords.clone();
    let center_in_original_space = basis * center_new_vec;
    let mut center = Point::new(center_in_original_space.iter().copied().collect());
    center.add(q0);

    Ball::new(center, miniball.radius)
}


fn extend_to_basis(points: &DMatrix<f64>) -> DMatrix<f64> {
    let dim = points.nrows();
    let n = points.ncols();

    let svd = points.clone().svd(true, false);
    let u = svd.u.unwrap();
    let sigma = svd.singular_values;

    let rank = sigma.iter().filter(|&&x| x > 1e-14).count();
    let basis = u.columns(0, rank).into_owned();

    basis
}


// generalised formula from https://mathworld.wolfram.com/Circumsphere.html
// https://en.wikipedia.org/wiki/Circumcircle
fn circumsphere<M> (boundary : &Vec<usize>, space : &PointCloud<M>) -> Ball
where 
    M: Euclidean
{
    let dim = space.dim();
    let n = dim + 1; // length of boundary

    let norms : Vec<f64> = (0..n).map(|x| space.norm_squared(space.get(boundary[x]))).collect();

    let mut c : Vec<f64> = Vec::new();
    for i in 0..dim {
        let M_i = DMatrix::from_fn(n, n, |row, col| {
                if col == 0 { return norms[row]; }
                if col == dim { return 1.0; }
                if col < i + 1 { space.get(boundary[row]).coords[col - 1] }
                else { space.get(boundary[row]).coords[col] }
        });
        c.push((-1.0_f64).powf(i as f64) * M_i.determinant());
    }

    let A = DMatrix::from_fn(n, n, |row, col| {
        if col == dim { return 1.0; }
        space.get(boundary[row]).coords[col]
    });
    let a = 1.0 / (2.0 * A.determinant());

    let mut center = Point::new(c);
    center.multiply(a);

    let radius = space.distance(&center, space.get(boundary[0]));

    Ball::new(center, radius)
}
