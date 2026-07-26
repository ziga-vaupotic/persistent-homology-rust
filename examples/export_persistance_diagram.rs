use persistent_homology::algebra::discretisation::{
    build_boundary_matrices, reduce_boundary_matrices,
};
use persistent_homology::algebra::persistence::compute_persistence_diagram;
use persistent_homology::construction::vietoris_rips;
use persistent_homology::geometry::EuclideanInnerProduct;
use persistent_homology::io::csv::{
    export_filtration_csv, export_persistence_csv, import_point_cloud_csv,
};

use std::path::Path;

fn main() {
    let path_name = "torus";
    let path_data = format!("examples/data/{}.csv", path_name);
    let path = Path::new(&path_data);

    let pointset = import_point_cloud_csv::<3, EuclideanInnerProduct>(path, EuclideanInnerProduct)
        .expect("Failed to read CSV");

    let filtration = vietoris_rips(&pointset, Some(4.0), Some(3));

    let path_filtration = format!("examples/data/{}_filtration.csv", path_name);
    export_filtration_csv(&path_filtration, &filtration).expect("Failed to export filtration!");

    let matrices = build_boundary_matrices(&filtration);

    let reduced_matrices = reduce_boundary_matrices(&matrices);

    let persistence = compute_persistence_diagram(&reduced_matrices);

    let path_persistence = format!("examples/data/{}_persistence.csv", path_name);

    export_persistence_csv(&path_persistence, &persistence).expect("Failed to export persistence!");
}
