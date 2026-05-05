use std::{error::Error, fs::File, path::Path};
use crate::geometry::point::Point;
use crate::geometry::point_set::PointSet;
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