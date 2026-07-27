use nalgebra::SVector;
use persistent_homology::construction::vietoris_rips;
use persistent_homology::geometry::{EuclideanInnerProduct, Point, PointCloud};

#[test]
fn test_filtration_is_sorted() {
    let points = vec![
        Point::<2>::new(SVector::<f64, 2>::from_row_slice(&[0.0, 0.0])),
        Point::<2>::new(SVector::<f64, 2>::from_row_slice(&[1.0, 0.0])),
        Point::<2>::new(SVector::<f64, 2>::from_row_slice(&[0.5, 1.0])),
    ];

    let pointset =
        PointCloud::new(points, EuclideanInnerProduct).expect("Pointset couldn't be generated.");

    let filtration = vietoris_rips(&pointset, None, Some(1));

    for i in 1..filtration.cells.len() {
        assert!(filtration.cells[i - 1].filtration_value <= filtration.cells[i].filtration_value);
    }
}

#[test]
fn test_simple_persistence_example() {
    // Two points close together, one far away
    let points = vec![
        Point::<2>::new(SVector::<f64, 2>::from_row_slice(&[0.0, 0.0])),
        Point::<2>::new(SVector::<f64, 2>::from_row_slice(&[0.1, 0.0])), // close to first
        Point::<2>::new(SVector::<f64, 2>::from_row_slice(&[10.0, 0.0])), // far away
    ];

    let pointset =
        PointCloud::new(points, EuclideanInnerProduct).expect("Pointset couldn't be generated.");

    let filtration = vietoris_rips(&pointset, None, Some(1));

    // Should have 3 points, 3 edges (0-1 close, 0-2 and 1-2 far)
    assert_eq!(filtration.cells.len(), 6); // 3 verts + 3 edges

    // The close edge should have lower filtration
    let close_edge_filtration = filtration
        .cells
        .iter()
        .find(|s| s.vertices == vec![0, 1])
        .unwrap()
        .filtration_value;
    let far_edge_filtration = filtration
        .cells
        .iter()
        .find(|s| s.vertices == vec![0, 2])
        .unwrap()
        .filtration_value;

    assert!(close_edge_filtration < far_edge_filtration);
}
