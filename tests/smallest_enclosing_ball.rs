use nalgebra::SVector;
use persistent_homology::geometry::{EuclideanInnerProduct, Point, PointCloud, seb::*};

use std::f64::consts::PI;

#[test]
fn in_2_dimensions() {
    let r = 10;
    let tolerance = 0.01;

    let mut point_set: Vec<Point<2>> = Vec::new();
    for _i in 0..r {
        // scatter random points inside the unit disk
        let (radius, phi) = (
            rand::random_range(0.0..=1.0),
            rand::random_range(0.0..(2.0 * PI)),
        );
        point_set.push(Point::<2>::new(SVector::<f64, 2>::from_row_slice(&[
            radius * phi.cos(),
            radius * phi.sin(),
        ])));
    }
    (0..2).for_each(|x| {
        point_set.push(Point::<2>::new(SVector::<f64, 2>::from_row_slice(&[
            (x as f64 * PI).cos(),
            (x as f64 * PI).sin(),
        ])))
    });

    let space = PointCloud::new(point_set, EuclideanInnerProduct).unwrap();

    let points: Vec<usize> = (0..(r + 2)).collect(); // every point

    let welzl_ball = welzl(&points, &space);
    let larsson_ball = larsson(&points, 0.01, &space);

    let center = Point::<2>::new(SVector::<f64, 2>::from_row_slice(&[0.0, 0.0]));
    let radius = 1.0;

    assert!((welzl_ball.o() - &center).is_zero());
    assert!((welzl_ball.r() - radius).abs() <= 1e-12);
    assert!(larsson_ball.r() <= radius * (1.0 + tolerance) + 1e-12);
}

#[test]
fn in_3_plus_dimensions() {
    let max_dim = 7;
    let tolerance = 0.01;

    for dim in 3..=max_dim {
        macro_rules! run_dim {
            ($d:literal) => {{
                let point_set: Vec<Point<$d>> =
                    (0..$d).map(|i| Point::<$d>::standard_unit(i)).collect();
                let space = PointCloud::new(point_set, EuclideanInnerProduct).unwrap();

                for i in 1..=$d {
                    let points: Vec<usize> = (0..i).collect();

                    let center = Point::<$d>::new(SVector::<f64, $d>::from_row_slice(
                        &(0..$d)
                            .map(|x| if x < i { 1.0 / i as f64 } else { 0.0 })
                            .collect::<Vec<f64>>(),
                    ));
                    let radius = space.distance(&center, space.get(0));

                    let welzl_ball = welzl(&points, &space);
                    let larsson_ball = larsson(&points, 0.01, &space);

                    assert!((welzl_ball.o() - &center).is_zero());
                    assert!((welzl_ball.r() - radius).abs() <= 1e-12);
                    assert!(larsson_ball.r() <= radius * (1.0 + tolerance));
                }
            }};
        }

        match dim {
            3 => run_dim!(3),
            4 => run_dim!(4),
            5 => run_dim!(5),
            6 => run_dim!(6),
            7 => run_dim!(7),
            _ => unreachable!(),
        }
    }
}
