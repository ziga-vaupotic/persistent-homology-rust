use plotters::prelude::*;
use vietoris_rips_rust::geometry::point_set::PointSet;
use vietoris_rips_rust::io::csv::{import_point_set, export_filtration_csv};
use vietoris_rips_rust::construction::rips::vietoris_rips;

use std::path::Path;


fn main() -> Result<(), Box<dyn std::error::Error>> {

    let path = std::path::Path::new("data.csv");

    let pointset = import_point_set::<2>(path)
        .expect("Failed to read CSV");

    let filtration = vietoris_rips(&pointset, 2);

    let root: DrawingArea<BitMapBackend<'_>, plotters::coord::Shift> = BitMapBackend::new("test.png", (640, 480)).into_drawing_area();
    root.fill(&WHITE);

    let root = root.margin(10, 10, 10, 10);
    // After this point, we should be able to construct a chart context
    let mut chart = ChartBuilder::on(&root)
        // Set the caption of the chart
        .caption("This is our first plot", ("sans-serif", 40).into_font())
        // Set the size of the label region
        .x_label_area_size(20)
        .y_label_area_size(40)
        // Finally attach a coordinate on the drawing area and make a chart context
        .build_cartesian_2d(-10f64..10f64, -10f64..10f64)?;



    chart.draw_series(
        PointSeries::of_element(
            (0..pointset.len()).map(|i| {
                let p = pointset.get(i);
                (p.coords[0], p.coords[1])
            }),
            2,
            &BLACK,
            &|c: (f64, f64), s, st| {
                EmptyElement::at(c)
                    + Circle::new((0, 0), s, st.filled())
            },
        )
    )?;
    root.present()?;
    Ok(())
}