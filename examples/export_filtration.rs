
use persistent_homology::algebra::discretisation::{ build_boundary_matrices, reduce_boundary_matrices };
use persistent_homology::algebra::persistence::compute_persistence_diagram;
use persistent_homology::geometry::{ EuclidianInnerProduct, PointCloud };
use persistent_homology::io::csv::{ export_filtration_csv, export_persistence_csv, import_point_set_csv };
use persistent_homology::construction::vietoris_rips;
use persistent_homology::construction::cech;


use std::path::Path;

fn main() {

    let path_name = "torus";
    let path_data = format!("examples/data/{}.csv", path_name);
    let path = std::path::Path::new(&path_data);

    let pointset = import_point_set_csv::<3, EuclidianInnerProduct>(path, EuclidianInnerProduct)
        .expect("Failed to read CSV");

    let filtration = cech(&pointset, Some(1.0), Some(2), 1e-6);

    let path_filtration = format!("examples/data/{}_filtration.csv", path_name);
    export_filtration_csv(&path_filtration, &filtration).expect("Failed to export filtration!");
}
