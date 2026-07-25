use plotters::prelude::*;
use persistent_homology::algebra::discretisation::{ build_boundary_matrices, reduce_boundary_matrices };
use persistent_homology::algebra::persistence::compute_persistence_diagram;
use persistent_homology::construction::vietoris_rips;
use persistent_homology::io::csv::import_point_set;

use std::path::Path;

struct PersistencePoint {
    birth: f64,
    death: f64,
    hom_dim: usize,
    infinite: bool,
}


fn main() -> Result<(), Box<dyn std::error::Error>> {
    let path = Path::new("data.csv");

    let pointset = import_point_set::<2>(path)?;
    let filtration = vietoris_rips(&pointset, None, Some(2));

    let matrices = build_boundary_matrices(&filtration);
    
    let mut diagram_points: Vec<PersistencePoint> = Vec::new();
    let mut max_filtration: f64 = 0.0;


    let matrices = build_boundary_matrices(&filtration);

    let reduced_matrices = reduce_boundary_matrices(&matrices);

    let persistence = compute_persistence_diagram(&reduced_matrices);


    for pair in &persistence.pairs {

        let (death, infinite) = match pair.death {
            Some(d) => (
                filtration.simplices[d].filtration_value,
                false,
            ),
            None => (
                f64::INFINITY, // or filtration.max_eps, etc.
                true,
            ),
        };

        diagram_points.push(PersistencePoint {
            birth: filtration.simplices[pair.birth].filtration_value,
            death: death,
            hom_dim: pair.dimension,
            infinite: infinite,
        });

    }

    let margin = 0.1 * max_filtration;
    let root: DrawingArea<BitMapBackend<'_>, plotters::coord::Shift> =
        BitMapBackend::new("persistence_diagram.png", (900, 900)).into_drawing_area();
    root.fill(&WHITE)?;

    let root = root.margin(50, 50, 50, 50);
    let mut chart = ChartBuilder::on(&root)
        .caption("Persistence Diagram", ("sans-serif", 32).into_font().style(FontStyle::Bold))
        .x_label_area_size(40)
        .y_label_area_size(50)
        .build_cartesian_2d(
            -margin..(max_filtration + margin),
            -margin..(max_filtration + margin),
        )?;

    chart
        .configure_mesh()
        .x_desc("Birth")
        .y_desc("Death")
        .draw()?;

    chart.draw_series(std::iter::once(PathElement::new(
        vec![(0.0, 0.0), (max_filtration, max_filtration)],
        ShapeStyle::from(&RGBColor(150, 150, 150)).stroke_width(2),
    )))?;

    let colors = vec![
        RGBColor(31, 119, 180),
        RGBColor(44, 160, 44),
        RGBColor(214, 39, 40),
        RGBColor(148, 103, 189),
        RGBColor(140, 86, 75),
    ];

    for dim in 0..matrices.len() {
        let points: Vec<_> = diagram_points
            .iter()
            .filter(|p| p.hom_dim == dim && !p.infinite)
            .collect();

        if points.is_empty() {
            continue;
        }

        let color = colors.get(dim).unwrap_or(&RGBColor(100, 100, 100));
        chart.draw_series(points.iter().map(|p| {
            Circle::new((p.birth, p.death), 1, color.filled())
        }))?;
    }

    for dim in 0..matrices.len() {
        let points: Vec<_> = diagram_points
            .iter()
            .filter(|p| p.hom_dim == dim && p.infinite)
            .collect();

        if points.is_empty() {
            continue;
        }

        let color = colors.get(dim).unwrap_or(&RGBColor(100, 100, 100));
        chart.draw_series(points.iter().map(|p| {
            TriangleMarker::new((p.birth, p.death), 1, color.filled())
        }))?;
    }

    root.present()?;
    println!("Persistence diagram saved to 'persistence_diagram.png'");
    for dim in 0..matrices.len() {
        let count = diagram_points.iter().filter(|p| p.hom_dim == dim).count();
        println!("H_{} features: {}", dim, count);
    }
    
    Ok(())
}
