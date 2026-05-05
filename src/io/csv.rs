use std::{error::Error, fs::File, path::Path, io::Write};
use crate::geometry::point::Point;
use crate::geometry::point_set::PointSet;

use crate::topology::filtration::Filtration;

use csv;

pub fn import_point_set<const D: usize>(
    path: &Path,
) -> Result<PointSet, Box<dyn Error>> {
    let file = File::open(path)?;
    let mut rdr = csv::ReaderBuilder::new()
        .has_headers(false)
        .from_reader(file);

    let mut points: Vec<Point> = Vec::new();

    for result in rdr.records() {
        let record = result?;

        if record.len() != D {
            return Err(format!("Expected {} columns, got {}", D, record.len()).into());
        }

        let mut arr = Vec::with_capacity(D);

        for v in record.iter() {
            arr.push(v.parse::<f64>()?);
        }

        points.push(Point { coords: arr });
    }

    Ok(PointSet::new(points)?)
}


pub fn export_filtration_csv(
    path: &str,
    filtration: &Filtration,
) -> Result<(), Box<dyn Error>> {
    let mut file = File::create(path)?;

    let mut simplices = filtration.simplices.clone();
    simplices.sort_by(|a, b| {
        a.filtration_value
            .partial_cmp(&b.filtration_value)
            .unwrap()
    });

    for simplex in simplices {
        write!(file, "{}", simplex.filtration_value)?;

        for v in simplex.vertices {
            write!(file, ",{}", v)?;
        }

        writeln!(file)?;
    }

    Ok(())
}