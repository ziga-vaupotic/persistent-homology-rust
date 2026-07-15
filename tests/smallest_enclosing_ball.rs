


use vietoris_rips_rust::geometry::{ Ball, Point, PointSet, seb::* };

use rand::prelude::SliceRandom;
use std::f64::consts::PI;


#[test]
fn in_2_dimensions() {
    let n = 3;
    let r = 10;
    let tolerance = 0.01;

    let mut point_set : Vec<Point> = Vec::new();
    for i in 0..n {
        let phi = 2.0 * PI * i as f64 / n as f64;
        point_set.push(Point::new(vec![phi.cos(), phi.sin()]));
    }
    for i in 0..r { // scatter random points inside the unit disk
        let (radius, phi) = (rand::random_range(0.0..=1.0), rand::random_range(0.0..(2.0 * PI)));
        point_set.push(Point::new(vec![radius * phi.cos(), radius * phi.sin()]));
    }
    point_set.shuffle(&mut rand::rng());
    let space = PointSet::new(point_set).unwrap();

    let center = Point::new(vec![0.0; 2]);
    let radius = 1.0;

    let points : Vec<usize> = (0..(n + r)).collect(); // every point

    let welzl_ball = welzl(&points, &space);
    let larsson_ball = larsson(&points, 0.01, &space);

    assert!(Point::difference(&welzl_ball.o(), &center).is_zero());
    assert!((welzl_ball.r() - radius).abs() <= 1e-12);
    assert!(larsson_ball.r() <= radius * (1.0 + tolerance) + 1e-12);
}


#[test]
fn in_3_plus_dimensions() {
    let max_dim = 7;
    let tolerance = 0.01;

    for dim in 3..=max_dim {
        let mut point_set : Vec<Point> = (0..dim).map(|i| Point::standard_unit(i, dim)).collect();
        let space = PointSet::new(point_set).unwrap();

        for i in 1..=dim {
            let points : Vec<usize> = (0..i).collect();

            let center = Point::new((0..dim).map(|x| if x < i { 1.0 / i as f64 } else { 0.0 }).collect());
            let radius = center.distance(space.get(0));

            let welzl_ball = welzl(&points, &space);
            let larsson_ball = larsson(&points, 0.01, &space);

            assert!(Point::difference(&welzl_ball.o(), &center).is_zero());
            assert!((welzl_ball.r() - radius).abs() <= 1e-12);
            assert!(larsson_ball.r() <= radius * (1.0 + tolerance));
        }
    }
}
