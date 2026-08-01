use nalgebra::SVector;
use persistent_homology::construction::{cech_exact, vietoris_rips};
use persistent_homology::geometry::{EuclideanCloud, Point, RealCoordinateSpace};

use std::f64::consts::PI;

#[test]
fn single_point() {
    let space = EuclideanCloud::new(vec![Point::new([0.0, 0.0, 0.0])]);

    let rips_filtration = vietoris_rips(&space, None, None);
    let cech_filtration = cech_exact(&space, None, None);

    assert!(rips_filtration.len() == 1);
    assert!(rips_filtration.cells[0].filtration_value == 0.0);

    assert!(cech_filtration.len() == 1);
    assert!(cech_filtration.cells[0].filtration_value == 0.0);
}

#[test]
fn triangle() {
    let n = 3;
    let mut point_set: Vec<Point<2>> = Vec::new();
    for i in 0..n {
        let phi = 2.0 * PI * i as f64 / n as f64;
        point_set.push(Point::new([phi.cos(), phi.sin()]));
    }
    let space = EuclideanCloud::new(point_set);

    let rips_filtration = vietoris_rips(&space, None, None);
    let cech_filtration = cech_exact(&space, None, None);

    assert!(rips_filtration.len() == 7);
    assert!(rips_filtration.max_dim() == 2);

    assert!(cech_filtration.len() == 7);
    assert!(cech_filtration.max_dim() == 2);

    for dim in 0..=2 {
        let rips_simplices = rips_filtration.cells_of_dim(dim);
        let cech_simplices = cech_filtration.cells_of_dim(dim);

        let i = dim as i32;
        let num = (-i.pow(3) + 2 * i.pow(2) - i + 3) as usize; // num(0) = 3, num(1) = 3, num(2) = 1

        assert!(rips_simplices.len() == num);
        assert!(cech_simplices.len() == num);

        let mut rips_value = 0.0;
        let mut cech_value = 0.0;
        match dim {
            1 => {
                rips_value = 2.0 * f64::cos(PI / 6.0);
                cech_value = f64::cos(PI / 6.0);
            }
            2 => {
                rips_value = 2.0 * f64::cos(PI / 6.0);
                cech_value = 1.0;
            }
            _ => {}
        }

        rips_simplices
            .iter()
            .for_each(|x| assert!((x.filtration_value - rips_value).abs() <= 1e-12));
        cech_simplices
            .iter()
            .for_each(|x| assert!((x.filtration_value - cech_value).abs() <= 1e-12));
    }
}

#[test]
fn grid_2_dim() {
    let n = 3;
    let mut point_set: Vec<Point<2>> = Vec::new();
    for i in 0..n {
        for j in 0..n {
            point_set.push(Point::new([i as f64, j as f64]));
        }
    }
    let space = EuclideanCloud::new(point_set);
    let rips_filtration = vietoris_rips(&space, Some(2.0_f64.sqrt()), None);
    let cech_filtration = cech_exact(&space, Some(2.0_f64.sqrt() / 2.0), None);

    assert!(rips_filtration.max_dim() == 3);
    assert!((rips_filtration.max_filtration_value() - 2.0_f64.sqrt()).abs() <= 1e-12);

    assert!(cech_filtration.max_dim() == 3);
    assert!((cech_filtration.max_filtration_value() - 2.0_f64.sqrt() / 2.0).abs() <= 1e-12);

    let dim_sizes = [9, 20, 16, 4];
    for (dim, x) in dim_sizes.iter().enumerate() {
        assert_eq!(rips_filtration.cells_of_dim(dim).len(), *x);
        assert_eq!(cech_filtration.cells_of_dim(dim).len(), *x);
    }
}

#[test]
fn grid_3_dim() {
    let n = 3;
    let mut point_set: Vec<Point<3>> = Vec::new();
    for i in 0..n {
        for j in 0..n {
            for k in 0..n {
                point_set.push(Point::new([i as f64, j as f64, k as f64]));
            }
        }
    }
    let space = EuclideanCloud::new(point_set);
    let rips_filtration = vietoris_rips(&space, Some(2.0_f64.sqrt()), None);
    let cech_filtration = cech_exact(&space, Some(2.0_f64.sqrt() / 2.0), None);

    assert!(rips_filtration.max_dim() == 3);
    assert!((rips_filtration.max_filtration_value() - 2.0_f64.sqrt()).abs() <= 1e-12);

    assert!(cech_filtration.max_dim() == 3);
    assert!((cech_filtration.max_filtration_value() - 2.0_f64.sqrt() / 2.0).abs() <= 1e-12);

    let rips_dim_sizes = [27, 126, 208, 116];
    let cech_dim_sizes = [27, 126, 144, 36];

    for dim in 0..=3 {
        assert_eq!(rips_filtration.cells_of_dim(dim).len(), rips_dim_sizes[dim]);
        assert_eq!(cech_filtration.cells_of_dim(dim).len(), cech_dim_sizes[dim]);
    }
}
