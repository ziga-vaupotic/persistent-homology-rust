use plotters::prelude::*;
use vietoris_rips_rust::io::csv::import_point_set;
use vietoris_rips_rust::construction::vietoris_rips;

use std::path::Path;


fn main() -> Result<(), Box<dyn std::error::Error>> {
    let path = Path::new("data.csv");

    let pointset = import_point_set::<2>(path)?;
    let filtration = vietoris_rips(&pointset, None, Some(2));
    let epsilon = 15.0;
    let complex_at_epsilon = filtration.get_simplicial_complex(epsilon);

    let xs: Vec<f64> = (0..pointset.len())
        .map(|i| pointset.get(i).coords[0])
        .collect();
    let ys: Vec<f64> = (0..pointset.len())
        .map(|i| pointset.get(i).coords[1])
        .collect();

    let x_min = xs.iter().cloned().fold(f64::INFINITY, f64::min);
    let x_max = xs.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let y_min = ys.iter().cloned().fold(f64::INFINITY, f64::min);
    let y_max = ys.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let margin = 1.0;

    let root: DrawingArea<BitMapBackend<'_>, plotters::coord::Shift> = BitMapBackend::new("test.png", (640, 480)).into_drawing_area();
    root.fill(&WHITE)?;

    let root = root.margin(10, 10, 10, 10);
    // After this point, we should be able to construct a chart context
    let mut chart = ChartBuilder::on(&root)
        .caption(
            format!("Vietoris–Rips complex at ε = {:.2}", epsilon),
            ("sans-serif", 28).into_font().style(FontStyle::Bold),
        )
        .x_label_area_size(30)
        .y_label_area_size(40)
        .build_cartesian_2d(
            (x_min - margin)..(x_max + margin),
            (y_min - margin)..(y_max + margin),
        )?;

    chart
        .configure_mesh()
        .x_desc("X")
        .y_desc("Y")
        .light_line_style(ShapeStyle::from(&RGBColor(200, 200, 200)).stroke_width(1))
        .draw()?;

    for triangle in complex_at_epsilon.simplices.iter().filter(|s| s.dimension() == 2) {
        let a = pointset.get(triangle.vertices[0]);
        let b = pointset.get(triangle.vertices[1]);
        let c = pointset.get(triangle.vertices[2]);
        chart.draw_series(std::iter::once(Polygon::new(
            vec![
                (a.coords[0], a.coords[1]),
                (b.coords[0], b.coords[1]),
                (c.coords[0], c.coords[1]),
            ],
            RGBColor(100, 170, 220).filled(),
        )))?;
    }

    for edge in complex_at_epsilon.simplices.iter().filter(|s| s.dimension() == 1) {
        let a = pointset.get(edge.vertices[0]);
        let b = pointset.get(edge.vertices[1]);
        chart.draw_series(LineSeries::new(
            vec![(a.coords[0], a.coords[1]), (b.coords[0], b.coords[1])],
            RGBColor(20, 90, 160).stroke_width(2),
        ))?;
    }

    chart.draw_series(
        PointSeries::of_element(
            (0..pointset.len()).map(|i| {
                let p = pointset.get(i);
                (p.coords[0], p.coords[1])
            }),
            4,
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
