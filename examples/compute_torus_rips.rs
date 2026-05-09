use vietoris_rips_rust::geometry::point_set::PointSet;
use vietoris_rips_rust::io::csv::{import_point_set, export_filtration_csv};
use vietoris_rips_rust::construction::rips::vietoris_rips;

use std::path::Path;

fn main() {
    let path = std::path::Path::new("data.csv");

    let pointset = import_point_set::<2>(path)
        .expect("Failed to read CSV");

    let filtration = vietoris_rips(&pointset, 2);

    export_filtration_csv("filtration.csv", &filtration).expect("Failed to export filtration!")
}