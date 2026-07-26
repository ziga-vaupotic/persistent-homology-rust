use crate::geometry::{Ball, Euclidean, Point, PointCloud};

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
pub fn larsson<const D: usize, M>(
    points: &[usize],
    epsilon: f64,
    space: &PointCloud<D, M>,
) -> Ball<D>
where
    M: Euclidean<D>,
{
    let delta = epsilon / 2.0;
    let mut C: Vec<Point<D>> = Vec::new();

    let (q_prime, _d) = farthest_point_space(space.get(points[0]), points, space);
    let (q, d) = farthest_point_space(&q_prime, points, space);

    let mut c = Ball::new(center_of_mass(&q, &q_prime), d / 2.0);

    C.push(q_prime);
    C.push(q);

    loop {
        let (p, h) = farthest_point_space(c.o(), points, space);
        if h <= c.r() * (1.0 + epsilon) {
            return c;
        }

        C.push(p.clone());
        c = update_ball(&c, &p, h);
        c = solve_approx_ball::<D, M>(c.clone(), &C, delta, space);
    }
}

fn solve_approx_ball<const D: usize, M>(
    mut c: Ball<D>,
    P: &[Point<D>],
    delta: f64,
    space: &PointCloud<D, M>,
) -> Ball<D>
where
    M: Euclidean<D>,
{
    loop {
        let (q, h) = farthest_point::<D, M>(c.o(), P, space);
        if h <= c.r() * (1.0 + delta) {
            return c;
        }
        c = update_ball(&c, &q, h);
    }
}

fn farthest_point<const D: usize, M>(
    p: &Point<D>,
    P: &[Point<D>],
    space: &PointCloud<D, M>,
) -> (Point<D>, f64)
where
    M: Euclidean<D>,
{
    let mut max_distance = 0.0;
    let mut farthest_point = p.clone();
    for x in P.iter() {
        let d = space.distance(p, x);
        if d >= max_distance {
            max_distance = d;
            farthest_point = x.clone();
        }
    }
    (farthest_point, max_distance)
}

fn farthest_point_space<const D: usize, M>(
    p: &Point<D>,
    P: &[usize],
    space: &PointCloud<D, M>,
) -> (Point<D>, f64)
where
    M: Euclidean<D>,
{
    let (mut max_distance, mut max_index) = (0.0, 0);
    for &i in P.iter() {
        let d = space.distance(space.get(i), p);
        if d >= max_distance {
            max_distance = d;
            max_index = i;
        }
    }
    (space.get(max_index).clone(), max_distance)
}

fn center_of_mass<const D: usize>(p: &Point<D>, q: &Point<D>) -> Point<D> {
    let mut o = p + q;
    o.multiply(1.0 / 2.0);
    o
}

fn update_ball<const D: usize>(old: &Ball<D>, p: &Point<D>, h: f64) -> Ball<D> {
    let mut a = old.o() - p;
    a.multiply(old.r() / h);
    Ball::new(&a + p, (old.r().powf(2.0) / h + h) / 2.0)
}
