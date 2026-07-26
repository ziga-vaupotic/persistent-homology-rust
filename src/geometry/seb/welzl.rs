use nalgebra::{DMatrix, SVector, DVector};
use rand::seq::SliceRandom;

use crate::geometry::{Ball, Euclidean, Point, PointCloud};

/// Compute the smallest enclosing ball using Welzl's algorithm (exact).
///
/// This is a randomized algorithm that exactly computes the smallest enclosing ball
/// of a set of points. It uses recursive subdivision and is guaranteed to find the optimal ball.
///
/// # Arguments
///
/// * `points` - Indices of points in the point cloud.
/// * `space` - The Euclidean point cloud.
///
/// # Returns
///
/// A `Ball` representing the smallest enclosing ball of the input points.
///
/// # Complexity
///
/// Expected $O((d+1)! \cdot n)$ where $d$ is dimension and $n$ is number of points.
/// Suitable for low dimensions (< 5, possibly up to 10).
///
/// # References
///
/// E. Welzl, Smallest enclosing disks (balls and ellipsoids), 1991
/// <https://doi.org/10.1007/BFb0038202>
///
/// # Note
///
/// For most practical purposes, consider using `larsson` for faster approximate computation.
pub fn welzl<const D: usize, M>(points: &[usize], space: &PointCloud<D, M>) -> Ball<D>
where
    M: Euclidean<D>,
{
    let mut P = points.to_owned();
    P.shuffle(&mut rand::rng());

    welzl_rec(&mut P, &mut Vec::new(), space)
}

// using Gärtner 1999 move-to-front implementation
// NOTE should probably be using LinkedList instead of Vec
// unfortunatelly rust LinkedLists do not support remove/insert operations (nightly versions only)
// TODO consider a move to VecDeque (for P only)
fn welzl_rec<const D: usize, M>(P: &mut [usize], B: &mut Vec<usize>, space: &PointCloud<D, M>) -> Ball<D>
where
    M: Euclidean<D>,
{
    let mut miniball = from_boundary(B, space);

    if P.is_empty() || B.len() == space.dim() + 1 {
        return miniball;
    }

    let n = P.len();
    for i in 0..n {
        let p = P[i];

        if space.contained_in_ball(&miniball, space.get(p)) {
            continue;
        }

        B.push(p);
        miniball = welzl_rec(&mut P[..i].to_vec(), B, space);
        B.pop();

        move_front(P, i);
    }
    miniball
}

fn move_front(v: &mut [usize], i: usize) {
    v[..=i].rotate_right(1);
}

fn from_boundary<const D: usize, M>(boundary: &[usize], space: &PointCloud<D, M>) -> Ball<D>
// |boundary| <= dim + 1
where
    M: Euclidean<D>,
{
    match boundary.len() {
        0 => Ball::new(Point::new(SVector::<f64, D>::zeros()), 0.0),
        1 => Ball::new(space.get(boundary[0]).clone(), 0.0),
        2 => {
            let mut o = space.get(boundary[0]) + space.get(boundary[1]);
            o.multiply(1.0 / 2.0);
            let r = space.distance(&o, space.get(boundary[0]));
            Ball::new(o, r)
        }
        x if x == space.dim() + 1 => {
            let points = boundary.iter().map(|&x| DVector::from_iterator(D, space.get(x).coords.iter().copied())).collect::<Vec<DVector<f64>>>();
            let (centre, radius) = circumsphere(&points);
            
            let centre_point = Point::new(
            SVector::<f64, D>::from_iterator(centre.iter().copied())
            );
            Ball::new(centre_point, radius)
        }
        _ => on_affine_subspace(boundary, space),
    }
}

// procedure :
// find smallest affine subspace containing boundary
// find isometry to R^n subset R^d, where n is the dimension of that subspace
// calculate miniball there then move the center back to original subspace
// no need to change the radius as we have an isometry
fn on_affine_subspace<const D: usize, M>(boundary: &[usize], space: &PointCloud<D, M>) -> Ball<D>
// 2 < |boundary| < dim + 1
where
    M: Euclidean<D>,
{
    let q0 = space.get(boundary[0]);
    let dim = space.dim();

    let linear_parts = DMatrix::from_fn(dim, boundary.len() - 1, |row, col| {
        space.get(boundary[col + 1]).coords[row] - q0.coords[row]
    });

    let basis = extend_to_basis(&linear_parts);
    let n = basis.ncols();

    let transformed = basis.transpose() * &linear_parts;

    let mut new_space_vectors: Vec<DVector<f64>> = (0..n)
        .map(|i| transformed.column(i).into_owned())
        .collect();

    new_space_vectors.push(DVector::zeros(transformed.nrows()));
    // find center in subpace spanned by {Q(q_i - q_0)}_i
    let (centre, radius) = circumsphere(&new_space_vectors);

    // center = basis * center_new + q_0
    let center_in_original_space = &basis * &centre;
    Ball::new(
        &Point::new(SVector::<f64, D>::from_iterator(center_in_original_space.iter().cloned()))
            + q0,
        radius,
    )
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
fn circumsphere(boundary: &Vec<DVector<f64>>) -> (DVector<f64>, f64)
{
    let dim = boundary.len() - 1;
    let n = dim + 1; // length of boundary

    let norms: Vec<f64> = (0..n)
        .map(|x| boundary[x].norm_squared())
        .collect();

    let mut c = DVector::<f64>::from_element(dim, 0.0);

    for i in 0..dim {
        // D_x, D_y, ... in wolfram reference
        let M_i = DMatrix::from_fn(n, n, |row, col| match col {
            0 => norms[row],
            x if x == dim => 1.0,
            _ => {
                if col < i + 1 {
                    return boundary[row][col - 1];
                }
                boundary[row][col]
            }
        });
        c[i] = (-1.0_f64).powi(i as i32) * M_i.determinant();
    }

    let A = DMatrix::from_fn(n, n, |row, col| {
        if col == dim {
            return 1.0;
        }
        boundary[row][col]
    });
    let a = 1.0 / (2.0 * A.determinant());

    let radius: f64 = (&c - &boundary[0]).norm();

   (c * a, radius)
}
