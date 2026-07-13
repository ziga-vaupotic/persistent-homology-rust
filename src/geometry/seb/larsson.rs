


use crate::geometry::{ Ball, Point, PointSet };


// T. Larsson, L. Källberg, Fast and Robust Approximation of Smallest Enclosing Balls
// in Arbitrary Dimensions, 2013
// https://doi.org/10.1111/cgf.12176
// NOTE O(dim n / epsilon + dim / epsilon^3) => good in higher dimensional spaces
pub fn larsson(points : &Vec<usize>, epsilon : f64, space : &PointSet) -> Ball {
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
        c = solve_approx_ball(c.clone(), &C, delta);
    }
}


fn solve_approx_ball(mut c : Ball, P : &Vec<Point>, delta : f64) -> Ball {
    loop {
        let (q, h) = farthest_point(c.o(), P);
        if h <= c.r() * (1.0 + delta) { return c; }
        c = update_ball(&c, &q, h);
    }
}


fn farthest_point(p : &Point, P : &Vec<Point>) -> (Point, f64) {
    let mut max_distance = 0.0;
    let mut farthest_point = p.clone();
    for x in P.iter() {
        let d = p.distance(x);
        if d >= max_distance {
            max_distance = d;
            farthest_point = x.clone();
        }
    }
    (farthest_point, max_distance)
}


fn farthest_point_space(p : &Point, P : &Vec<usize>, space : &PointSet) -> (Point, f64) {
    let (mut max_distance, mut max_index) = (0.0, 0);
    for &i in P.iter() {
        let d = space.get(i).distance(p);
        if d >= max_distance {
            max_distance = d;
            max_index = i;
        }
    }
    (space.get(max_index).clone(), max_distance)
}


fn center_of_mass(p : &Point, q : &Point) -> Point {
    let mut o = p.clone();
    o.add(q);
    o.multiply(1.0 / 2.0);
    o
}


fn update_ball(old : &Ball, p : &Point, h : f64) -> Ball {
    let mut new = Point::difference(old.o(), p);
    new.multiply(old.r() / h);
    new.add(p);
    Ball::new(new, (old.r().powf(2.0) / h + h) / 2.0)
}
