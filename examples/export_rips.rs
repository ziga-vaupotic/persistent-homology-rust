


use vietoris_rips_rust::algebra::discretisation::{ build_boundary_matrices, reduce_boundary_matrices };
use vietoris_rips_rust::algebra::persistence::compute_persistence_diagram;
use vietoris_rips_rust::geometry::PointSet;
use vietoris_rips_rust::io::csv::{ export_filtration_csv, export_persistence_csv, import_point_set };
use vietoris_rips_rust::construction::vietoris_rips;
use vietoris_rips_rust::construction::cech;


use std::path::Path;

fn main() {

    let path_name = "circle";
    let path_data = format!("examples/data/{}.csv", path_name);
    let path = std::path::Path::new(&path_data);

    let pointset = import_point_set::<2>(path)
        .expect("Failed to read CSV");

    let filtration = cech(&pointset, None, Some(2), 1e-6);



    let path_filtration = format!("examples/data/{}_filtration.csv", path_name);
    export_filtration_csv(&path_filtration, &filtration).expect("Failed to export filtration!");

    let matrices = build_boundary_matrices(&filtration);

    let reduced_matrices = reduce_boundary_matrices(&matrices);

    let persistence = compute_persistence_diagram(&reduced_matrices);

    let path_persistence = format!("examples/data/{}_persistence.csv", path_name);

    export_persistence_csv(&path_persistence, &persistence).expect("Failed to export persistence!");
}
