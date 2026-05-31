use vietoris_rips_rust::algebra::discretisation::{build_boundary_matrices, reduce_boundary_matrices};
use vietoris_rips_rust::algebra::persistence::{compute_persistence_diagram};
use vietoris_rips_rust::geometry::point_set::PointSet;
use vietoris_rips_rust::io::csv::{export_filtration_csv, export_persistence_csv, import_point_set};
use vietoris_rips_rust::construction::rips::vietoris_rips;

use std::path::Path;

fn main() {
    let path = std::path::Path::new("data.csv");

    let pointset = import_point_set::<2>(path)
        .expect("Failed to read CSV");

    let filtration = vietoris_rips(&pointset, 2, None);

    export_filtration_csv("filtration.csv", &filtration).expect("Failed to export filtration!");

    let matrices = build_boundary_matrices(&filtration);

    let reduced_matrices = reduce_boundary_matrices(&matrices);

    let persistence = compute_persistence_diagram(&reduced_matrices);

    export_persistence_csv("persistence.csv", &persistence).expect("Failed to export persistence!");
}