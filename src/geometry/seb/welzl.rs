use nalgebra::DMatrix;
use rand::seq::SliceRandom;

use crate::geometry::{Ball, Euclidean, Point, PointCloud};

// E. Welzl, Smallest enclosing disks (balls and ellipsoids), 1991
// https://doi.org/10.1007/BFb0038202
// NOTE should be only used for low dimensions ie. < 5 (maybe up to 10? with a very small space)
// the basic algorith listed in the article is O((dim + 1)(dim + 1)!n), while time complexity for
// this specific version is not provided, it cannot be much better
// this implementation also calculates the ball from a given boundary from scratch every
// time so Gärtner 1999 should be prefered in most cases (once its implemented)
pub fn welzl<M>(points: &[usize], space: &PointCloud<M>) -> Ball
where
    M: Euclidean,
{
    let mut P = points.to_owned();
    P.shuffle(&mut rand::rng());

    welzl_rec(&mut P, &mut Vec::new(), space)
}

// using Gärtner 1999 move-to-front implementation
// NOTE should probably be using LinkedList instead of Vec
// unfortunatelly rust LinkedLists do not support remove/insert operations (nightly versions only)
// TODO consider a move to VecDeque (for P only)
fn welzl_rec<M>(P: &mut Vec<usize>, B: &mut Vec<usize>, space: &PointCloud<M>) -> Ball
where
    M: Euclidean,
{
    let mut miniball = from_boundary(&B, space);
    if P.is_empty() || B.len() == space.dim() + 1 {
        return miniball;
    }

    let n = P.len();
    for i in 0..n {
        if space.contained_in_ball(&miniball, space.get(P[i])) {
            continue;
        }

        B.push(P[i]);
        miniball = welzl_rec(&mut P[..i].to_vec(), B, space);
        B.pop();

        move_front(P, i);
    }
    miniball
}

fn move_front(v: &mut Vec<usize>, x: usize) {
    let u = v[x];
    v.remove(x);
    v.splice(..0, [u]);
}

fn from_boundary<M>(boundary: &[usize], space: &PointCloud<M>) -> Ball
// |boundary| <= dim + 1
where
    M: Euclidean,
{
    match boundary.len() {
        0 => Ball::new(Point::new(Vec::new()), 0.0),
        1 => Ball::new(space.get(boundary[0]).clone(), 0.0),
        2 => {
            let mut o = space.get(boundary[0]) + space.get(boundary[1]);
            o.multiply(1.0 / 2.0);
            let r = space.distance(&o, space.get(boundary[0]));
            Ball::new(o, r)
        }
        x if x == space.dim() + 1 => circumsphere(boundary, space),
        _ => on_affine_subspace(boundary, space),
    }
}

// procedure :
// find smallest affine subspace containing boundary
// find isometry to R^n subset R^d, where n is the dimension of that subspace
// calculate miniball there then move the center back to original subspace
// no need to change the radius as we have an isometry
fn on_affine_subspace<M>(boundary: &[usize], space: &PointCloud<M>) -> Ball
// 2 < |boundary| < dim + 1
where
    M: Euclidean,
{
    let q0 = space.get(boundary[0]);
    let dim = space.dim();

    let linear_parts = DMatrix::from_fn(dim, boundary.len() - 1, |row, col| {
        space.get(boundary[col + 1]).coords[row] - q0.coords[row]
    });

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
    let new_boundary: Vec<usize> = (0..(n + 1)).collect();

    // find center in subpace spanned by {Q(q_i - q_0)}_i
    let miniball = circumsphere(&new_boundary, &new_space);

    // center = basis * center_new + q_0
    let center_new_vec = miniball.o().coords.clone();
    let center_in_original_space = basis * center_new_vec;
    Ball::new(&Point::new(center_in_original_space) + q0, miniball.radius)
}

fn extend_to_basis(points: &DMatrix<f64>) -> DMatrix<f64> {
    let svd = points.clone().svd(true, false);
    let u = svd.u.unwrap();
    let sigma = svd.singular_values;

    let rank = sigma.iter().filter(|&&x| x > 1e-14).count();

    u.columns(0, rank).into_owned()
}

// generalised formula from https://mathworld.wolfram.com/Circumsphere.html
// https://en.wikipedia.org/wiki/Circumcircle
fn circumsphere<M>(boundary: &[usize], space: &PointCloud<M>) -> Ball
where
    M: Euclidean,
{
    let dim = space.dim();
    let n = dim + 1; // length of boundary

    let norms: Vec<f64> = (0..n)
        .map(|x| space.norm_squared(space.get(boundary[x])))
        .collect();

    let mut c: Vec<f64> = Vec::new();
    for i in 0..dim {
        // D_x, D_y, ... in wolfram reference
        let M_i = DMatrix::from_fn(n, n, |row, col| match col {
            0 => norms[row],
            x if x == dim => 1.0,
            _ => {
                if col < i + 1 {
                    return space.get(boundary[row]).coords[col - 1];
                }
                space.get(boundary[row]).coords[col]
            }
        });
        c.push((-1.0_f64).powf(i as f64) * M_i.determinant());
    }

    let A = DMatrix::from_fn(n, n, |row, col| {
        if col == dim {
            return 1.0;
        }
        space.get(boundary[row]).coords[col]
    });
    let a = 1.0 / (2.0 * A.determinant());

    let mut center = Point::new(c);
    center.multiply(a);

    let radius = space.distance(&center, space.get(boundary[0]));

    Ball::new(center, radius)
}
