use crate::geometry::{Point, PointCloud, Set, Space};
use nalgebra::SVector;
use std::{error::Error, fs::File, io::Write, path::Path};

use crate::algebra::persistence::PersistenceDiagram;
use crate::topology::{Filtration, Simplex};

use csv;

/// Import a point cloud from a CSV file.
///
/// Reads a CSV file where each row is a point and each column is a coordinate.
/// The generic parameter `D` specifies the expected coordinate dimension of each point.
///
/// # Arguments
///
/// * `path` - Path to the CSV file.
/// * `geometry` - The metric or inner product space for the point cloud.
///
/// # Returns
///
/// A `PointCloud` if successful, or an error if the file cannot be read or
/// if points have inconsistent dimensions.
///
/// # Example
///
/// ```ignore
/// use persistent_homology::io::csv::import_point_cloud_csv;
/// use persistent_homology::geometry::EuclideanInnerProduct;
/// use std::path::Path;
///
/// let cloud = import_point_cloud_csv::<2, _>(Path::new("points.csv"), EuclideanInnerProduct)?;
/// ```
pub fn import_point_cloud_csv<const N: usize, S: Space>(
    path: &Path,
) -> Result<PointCloud<S>, Box<dyn Error>>
where
    S::Set: Set<Element = Point<N>>,
{
    let file = File::open(path)?;
    let mut rdr = csv::ReaderBuilder::new()
        .has_headers(false)
        .from_reader(file);

    let mut points: Vec<Point<N>> = Vec::new();

    for result in rdr.records() {
        let record = result?;

        if record.len() != N {
            return Err(format!("Expected {} columns, got {}", N, record.len()).into());
        }

        let mut arr = Vec::with_capacity(N);

        for v in record.iter() {
            arr.push(v.parse::<f64>()?);
        }

        points.push(Point::new(SVector::<f64, N>::from_row_slice(&arr)));
    }

    Ok(PointCloud::new(points))
}

/// Export a filtration to a CSV file.
///
/// Writes each simplex to a row with columns for filtration value, dimension, and vertices.
/// Simplices are sorted by filtration value.
///
/// # Arguments
///
/// * `path` - Path where the CSV file will be written.
/// * `filtration` - The filtration to export.
///
/// # Returns
///
/// `Ok(())` if successful, or an error if the file cannot be written.
pub fn export_filtration_csv(
    path: &str,
    filtration: &Filtration<Simplex>,
) -> Result<(), Box<dyn Error>> {
    let mut file = File::create(path)?;

    let simplices = filtration.cells.clone();

    for simplex in simplices {
        write!(file, "{}", simplex.filtration_value)?;

        for v in simplex.vertices {
            write!(file, ",{}", v)?;
        }

        writeln!(file)?;
    }

    Ok(())
}

/// Export a persistence diagram to a CSV file.
///
/// Writes dimension of homological feature, birth and death
///
/// # Arguments
///
/// * `path` - Path where the CSV file will be written.
/// * `persistence` - The Persistence diagram to export.
///
/// # Returns
///
/// `Ok(())` if successful, or an error if the file cannot be written.
pub fn export_persistence_csv(
    path: &str,
    persistence: &PersistenceDiagram,
) -> Result<(), Box<dyn Error>> {
    let mut file = File::create(path)?;

    for pair in &persistence.pairs {
        write!(file, "{},", pair.dimension)?;
        write!(file, "{},", pair.birth)?;
        match pair.death {
            Some(death) => write!(file, "{}", death)?,
            None => write!(file, "")?, // Infinite death: leave empty
        }
        writeln!(file)?;
    }
    Ok(())
}
