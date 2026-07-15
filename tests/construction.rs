


use vietoris_rips_rust::geometry::{ Point, PointSet };
use vietoris_rips_rust::construction::{ cech, cech_exact, vietoris_rips };

use std::f64::consts::PI;


#[test]
fn single_point() {
    let space = PointSet::new(vec![Point::new(vec![0.0; 3])]).unwrap();

    let rips_filtration = vietoris_rips(&space, None, None);
    let cech_filtration = cech_exact(&space, None, None);

    assert!(rips_filtration.size() == 1 && cech_filtration.size() == 1);
    assert!(rips_filtration.simplices[0].filtration_value == 0.0);
    assert!(cech_filtration.simplices[0].filtration_value == 0.0);
}


#[test]
fn triangle() {
    let n = 3;
    let mut point_set : Vec<Point> = Vec::new();
    for i in 0..n {
        let phi = 2.0 * PI * i as f64 / n as f64;
        point_set.push(Point::new(vec![phi.cos(), phi.sin()]));
    }
    let space = PointSet::new(point_set).unwrap();

    let rips_filtration = vietoris_rips(&space, None, None);
    let cech_filtration = cech_exact(&space, None, None);

    assert!(rips_filtration.size() == 7 && cech_filtration.size() == 7);
    assert!(rips_filtration.max_dim() == 2 && cech_filtration.max_dim() == 2);

    for dim in 0..=2 {
        let rips_simplices = rips_filtration.simplices_of_dim(dim);
        let cech_simplices = cech_filtration.simplices_of_dim(dim);
        
        let i = dim as i32;
        let num = (-i.pow(3) + 2 * i.pow(2) - i + 3) as usize; // num(0) = 3, num(1) = 3, num(2) = 1

        assert!(rips_simplices.len() == num && cech_simplices.len() == num);

        let mut rips_value = 0.0;
        let mut cech_value = 0.0;
        match dim {
            1 => { rips_value = 2.0 * f64::cos(PI / 6.0); cech_value = f64::cos(PI / 6.0); },
            2 => { rips_value = 2.0 * f64::cos(PI / 6.0); cech_value = 1.0; },
            _ => {}
        }

        for x in rips_simplices { assert!((x.filtration_value - rips_value).abs() <= 1e-12); }
        for x in cech_simplices { assert!((x.filtration_value - cech_value).abs() <= 1e-12); }
    }
}


#[test]
fn grid_2_dim() {
    let n = 3;
    let mut point_set : Vec<Point> = Vec::new();
    for i in 0..n {
        for j in 0..n {
            point_set.push(Point::new(vec![i as f64, j as f64]));
        }
    }
    let space = PointSet::new(point_set).unwrap();
    let rips_filtration = vietoris_rips(&space, Some(2.0_f64.sqrt()), None);
    let cech_filtration = cech_exact(&space, Some(2.0_f64.sqrt() / 2.0), None);

    assert!(rips_filtration.max_dim() == 3 && cech_filtration.max_dim() == 3);
    assert!((rips_filtration.max_filtration_value() - 2.0_f64.sqrt()).abs() <= 1e-12);
    assert!((cech_filtration.max_filtration_value() - 2.0_f64.sqrt() / 2.0).abs() <= 1e-12);

    let dim_sizes = [9, 20, 16, 4];
    for dim in 0..=3 {
        assert_eq!(rips_filtration.simplices_of_dim(dim).len(), dim_sizes[dim]);
        assert_eq!(cech_filtration.simplices_of_dim(dim).len(), dim_sizes[dim]);
    }
}


#[test]
fn grid_3_dim() {
    let n = 3;
    let mut point_set : Vec<Point> = Vec::new();
    for i in 0..n {
        for j in 0..n {
            for k in 0..n {
                point_set.push(Point::new(vec![i as f64, j as f64, k as f64]));
            }
        }
    }
    let space = PointSet::new(point_set).unwrap();
    let rips_filtration = vietoris_rips(&space, Some(2.0_f64.sqrt()), None);
    let cech_filtration = cech_exact(&space, Some(2.0_f64.sqrt() / 2.0), None);

    assert!(rips_filtration.max_dim() == 3 && cech_filtration.max_dim() == 3);
    assert!((rips_filtration.max_filtration_value() - 2.0_f64.sqrt()).abs() <= 1e-12);
    assert!((cech_filtration.max_filtration_value() - 2.0_f64.sqrt() / 2.0).abs() <= 1e-12);

    let rips_dim_sizes = [27, 126, 208, 116];
    let cech_dim_sizes = [27, 126, 144, 36];

    for dim in 0..=3 {
        assert_eq!(rips_filtration.simplices_of_dim(dim).len(), rips_dim_sizes[dim]);
        assert_eq!(cech_filtration.simplices_of_dim(dim).len(), cech_dim_sizes[dim]);
    }
}
