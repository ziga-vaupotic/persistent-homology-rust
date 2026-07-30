use persistent_homology::construction::cech;
use persistent_homology::io::csv::{export_filtration_csv, import_point_cloud_csv};

use std::path::Path;

fn main() {
    let path_name = "figure_eight";
    let path_data = format!("examples/data/{}.csv", path_name);
    let path = Path::new(&path_data);

    let pointset = import_point_cloud_csv::<2>(path)
        .expect("Failed to read CSV");

    let filtration = cech(&pointset, Some(0.5), Some(2), 1e-6);

    let path_filtration = format!("examples/data/{}_filtration.csv", path_name);
    export_filtration_csv(&path_filtration, &filtration).expect("Failed to export filtration!");
}
