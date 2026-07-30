use crate::geometry::{Ball, EuclideanBall, EuclideanCloud, EuclideanSpace, MetricSpace, Point};

/// Compute an approximate smallest enclosing ball using Larsson's algorithm.
///
/// This is a fast approximate algorithm for computing the smallest enclosing ball
/// in arbitrary dimensions. It provides a guaranteed approximation ratio of $(1 + \epsilon)$.
///
/// # Arguments
///
/// * `points` - Indices of points in the point cloud.
/// * `epsilon` - Approximation tolerance. The algorithm guarantees:
///   $$r_{\text{optimal}} \leq r_{\text{larsson}} \leq r_{\text{optimal}} (1 + \epsilon)$$
/// * `space` - The Euclidean point cloud.
///
/// # Returns
///
/// A `Ball` that approximately encloses all input points.
///
/// # Complexity
///
/// $O(d \cdot n / \epsilon + d / \epsilon^3)$ where $d$ is dimension and $n$ is number of points.
/// Performs well in higher-dimensional spaces.
///
/// # References
///
/// T. Larsson, L. Källberg, Fast and Robust Approximation of Smallest Enclosing Balls
/// in Arbitrary Dimensions, 2013
/// <https://doi.org/10.1111/cgf.12176>
pub fn larsson<const N: usize>(
    points: &[usize],
    epsilon: f64,
    space: &EuclideanCloud<N>,
) -> EuclideanBall<N> {
    let delta = epsilon / 2.0;
    let mut C: Vec<Point<N>> = Vec::new();

    let (q_prime, _d) = farthest_point_space(space.get(points[0]), points, space);
    let (q, d) = farthest_point_space(&q_prime, points, space);

    let mut c = Ball::new(center_of_mass(&q, &q_prime), d / 2.0);

    C.push(q_prime);
    C.push(q);

    loop {
        let (p, h) = farthest_point_space(c.o(), points, space);
        if h <= c.r() * (1.0 + epsilon) {
            return c;
            return EuclideanBall::new(c.center, h);
        }

        C.push(p.clone());
        c = update_ball(&c, &p, h);
        c = solve_approx_ball(c.clone(), &C, delta);
    }
}

fn solve_approx_ball<const N: usize>(
    mut c: Ball<EuclideanSpace<N>>,
    P: &[Point<N>],
    delta: f64,
) -> EuclideanBall<N> {
    loop {
        let (q, h) = farthest_point(c.o(), P);
        if h <= c.r() * (1.0 + delta) {
            return c;
        }
        c = update_ball(&c, &q, h);
    }
}

fn farthest_point<const N: usize>(p: &Point<N>, P: &[Point<N>]) -> (Point<N>, f64) {
    let mut max_distance = 0.0;
    let mut farthest_point = p.clone();
    for x in P.iter() {
        let d = EuclideanSpace::distance(p, x);
        if d >= max_distance {
            max_distance = d;
            farthest_point = x.clone();
        }
    }
    (farthest_point, max_distance)
}

fn farthest_point_space<const N: usize>(
    p: &Point<N>,
    P: &[usize],
    space: &EuclideanCloud<N>,
) -> (Point<N>, f64) {
    let (mut max_distance, mut max_index) = (0.0, 0);
    for &i in P.iter() {
        let d = EuclideanSpace::distance(space.get(i), p);
        if d >= max_distance {
            max_distance = d;
            max_index = i;
        }
    }
    (space.get(max_index).clone(), max_distance)
}

fn center_of_mass<const D: usize>(p: &Point<D>, q: &Point<D>) -> Point<D> {
    (1.0 / 2.0) * &(p + q)
}

fn update_ball<const N: usize>(
    old: &Ball<EuclideanSpace<N>>,
    p: &Point<N>,
    h: f64,
) -> EuclideanBall<N> {
    Ball::new(
        &((old.r() / h) * &(old.o() - p)) + p,
        (old.r().powf(2.0) / h + h) / 2.0,
    )
}
