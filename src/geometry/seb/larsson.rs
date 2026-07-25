


use crate::geometry::{ Ball, Point, PointCloud, Euclidean };


// T. Larsson, L. Källberg, Fast and Robust Approximation of Smallest Enclosing Balls
// in Arbitrary Dimensions, 2013
// https://doi.org/10.1111/cgf.12176
// NOTE O(dim n / epsilon + dim / epsilon^3) => good in higher dimensional spaces
// NOTE r_optimal <= r_larsson <= r_optimal (1 + epsilon)
pub fn larsson<M>(points : &Vec<usize>, epsilon : f64, space : &PointCloud<M>) -> Ball
where
    M : Euclidean
{
    let delta = epsilon / 2.0;
    let mut C : Vec<Point> = Vec::new();

    let (q_prime, _d) = farthest_point_space(space.get(points[0]), points, space);
    let (q, d) = farthest_point_space(&q_prime, points, space);

    let mut c = Ball::new(center_of_mass(&q, &q_prime), d / 2.0);

    C.push(q_prime);
    C.push(q);

    loop {
        let (p, h) = farthest_point_space(c.o(), points, space);
        if h <= c.r() * (1.0 + epsilon) { return c; }

        C.push(p.clone());
        c = update_ball(&c, &p, h);
        c = solve_approx_ball::<M>(c.clone(), &C, delta, space);
    }
}


fn solve_approx_ball<M>(mut c : Ball, P : &Vec<Point>, delta : f64, space : &PointCloud<M>) -> Ball
where
    M : Euclidean
{
    loop {
        let (q, h) = farthest_point::<M>(c.o(), P, space);
        if h <= c.r() * (1.0 + delta) { return c; }
        c = update_ball(&c, &q, h);
    }
}


fn farthest_point<M>(p : &Point, P : &Vec<Point>, space : &PointCloud<M>) -> (Point, f64)
where
    M : Euclidean
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


fn farthest_point_space<M>(p : &Point, P : &Vec<usize>, space : &PointCloud<M>) -> (Point, f64)
where
    M : Euclidean
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


fn center_of_mass(p : &Point, q : &Point) -> Point {
    let mut o = p + q;
    o.multiply(1.0 / 2.0);
    o
}


fn update_ball(old : &Ball, p : &Point, h : f64) -> Ball {
    let mut a = old.o() - p;
    a.multiply(old.r() / h);
    Ball::new(&a + p, (old.r().powf(2.0) / h + h) / 2.0)
}
