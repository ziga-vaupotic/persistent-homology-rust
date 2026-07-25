use persistent_homology::geometry::{ EuclidianInnerProduct, Point, PointCloud };
use persistent_homology::construction::vietoris_rips;


#[test]
fn test_filtration_is_sorted() {
    let points = vec![
        Point::new(vec![0.0, 0.0]),
        Point::new(vec![1.0, 0.0]),
        Point::new(vec![0.5, 1.0]),
    ];

    let pointset = PointCloud::new(points, EuclidianInnerProduct).expect("Pointset couldn't be generated.");

    let filtration = vietoris_rips(&pointset, None, Some(1));

    for i in 1..filtration.simplices.len() {
        assert!(filtration.simplices[i - 1].filtration_value <= filtration.simplices[i].filtration_value);
    }
}


#[test]
fn test_simple_persistence_example() {
    // Two points close together, one far away
    let points = vec![
        Point::new(vec![0.0, 0.0]),
        Point::new(vec![0.1, 0.0]), // close to first
        Point::new(vec![10.0, 0.0]), // far away
    ];

    let pointset = PointCloud::new(points, EuclidianInnerProduct).expect("Pointset couldn't be generated.");

    let filtration = vietoris_rips(&pointset, None, Some(1));

    // Should have 3 points, 3 edges (0-1 close, 0-2 and 1-2 far)
    assert_eq!(filtration.simplices.len(), 6); // 3 verts + 3 edges

    // The close edge should have lower filtration
    let close_edge_filtration = filtration.simplices.iter()
        .find(|s| s.vertices == vec![0, 1])
        .unwrap()
        .filtration_value;
    let far_edge_filtration = filtration.simplices.iter()
        .find(|s| s.vertices == vec![0, 2])
        .unwrap()
        .filtration_value;

    assert!(close_edge_filtration < far_edge_filtration);
}
