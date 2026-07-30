use crate::construction::{Construction, cliques};
use crate::geometry::{MetricSpace, PointCloud};
use crate::topology::{Filtration, Simplex};

use itertools::Itertools;

/// Construct a Vietoris-Rips complex from a point cloud.
///
/// The Vietoris-Rips complex includes a $k$-simplex if all pairwise distances among its vertices
/// are at most `epsilon`. The filtration value of each simplex is the maximum pairwise distance
/// in that simplex.
///
/// # Arguments
///
/// * `space` - The point cloud with an associated metric.
/// * `max_epsilon` - Maximum filtration value. Defaults to infinity if `None`.
/// * `max_dim` - Maximum simplex dimension. Defaults to no limit if `None`.
///
/// # Returns
///
/// A `Filtration` ordered by filtration value, containing all Vietoris-Rips simplices.
///
/// # Example
///
/// ```ignore
/// use persistent_homology::construction::vietoris_rips;
/// use nalgebra::SVector;
/// use persistent_homology::geometry::{Point, PointCloud, EuclideanInnerProduct};
///
/// let points = vec![
///     Point::<1>::new(SVector::<f64, 1>::from_row_slice(&[0.0])),
///     Point::<1>::new(SVector::<f64, 1>::from_row_slice(&[1.0])),
/// ];
/// let cloud = PointCloud::new(points, EuclideanInnerProduct).unwrap();
/// let filtration = vietoris_rips(&cloud, Some(2.0), Some(2));
/// ```
pub fn vietoris_rips<M: MetricSpace>(
    cloud: &PointCloud<M>,
    max_epsilon: Option<f64>,
    max_dim: Option<usize>,
) -> Filtration<Simplex>
{
    let max_epsilon = max_epsilon.unwrap_or(f64::MAX);
    let max_dim = max_dim.unwrap_or(usize::MAX - 1);

    let mut cons = Construction::new(max_dim, max_epsilon, 0.0, cloud);
    if max_dim == 0 {
        return Filtration::new(cons.simplices);
    }

    if !cons.traverse_edges(cloud, 1.0, true) {
        return Filtration::new(cons.simplices);
    }

    cliques::find_all(cloud, rips_radius, &mut cons);

    Filtration::new(cons.simplices)
}

fn rips_radius<M: MetricSpace>(
    clique: &[usize],
    _space: &PointCloud<M>,
    cons: &Construction,
) -> Option<f64>
{
    let mut max_d = 0.0;
    for v in (0..clique.len()).combinations(2) {
        let d = cons.distance[&(clique[v[0]], clique[v[1]])];
        max_d = if d >= max_d { d } else { max_d };
    }
    Some(max_d)
}
