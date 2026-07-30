use crate::construction::{Construction, cliques};
use crate::geometry::{EuclideanSpace, EuclideanCloud, seb};
use crate::topology::{Filtration, Simplex};

/// Construct a Čech complex from a point cloud.
///
/// The Čech complex includes a $k$-simplex if the smallest enclosing ball (SEB) of its vertices
/// has radius at most `epsilon`. The filtration value is the radius of the SEB.
///
/// # Arguments
///
/// * `space` - The point cloud with a Euclidean inner product space.
/// * `max_epsilon` - Maximum filtration value. Defaults to infinity if `None`.
/// * `max_dim` - Maximum simplex dimension. Defaults to no limit if `None`.
/// * `radius_tolerance` - Tolerance for SEB computation. Set to 0 for exact computation.
///   Use a positive value for faster approximate computation with the Larsson algorithm.
///
/// # Returns
///
/// A `Filtration` ordered by filtration value, containing all Čech simplices.
///
/// # Example
///
/// ```ignore
/// use persistent_homology::construction::cech;
/// use nalgebra::SVector;
/// use persistent_homology::geometry::{Point, PointCloud, EuclideanInnerProduct};
///
/// let points = vec![
///     Point::<2>::new(SVector::<f64, 2>::from_row_slice(&[0.0, 0.0])),
///     Point::<2>::new(SVector::<f64, 2>::from_row_slice(&[1.0, 0.0])),
/// ];
/// let cloud = PointCloud::new(points, EuclideanInnerProduct).unwrap();
/// let filtration = cech(&cloud, Some(2.0), Some(2), 0.0);
/// ```
pub fn cech<const N: usize>(
    cloud: &EuclideanCloud<N>,
    max_epsilon: Option<f64>,
    max_dim: Option<usize>,
    radius_tolerance: f64,
) -> Filtration<Simplex>
{
    let max_epsilon = max_epsilon.unwrap_or(f64::MAX / 2.0);
    let max_dim = max_dim.unwrap_or(usize::MAX - 1);
    let radius_tolerance = radius_tolerance.abs();

    let mut cons = Construction::new(max_dim, max_epsilon, radius_tolerance, cloud);
    if max_dim == 0 {
        return Filtration::new(cons.simplices);
    }

    if !cons.traverse_edges(cloud, 2.0, false) {
        return Filtration::new(cons.simplices);
    }

    let in_ball = if radius_tolerance == 0.0 {
        in_ball_exact
    } else {
        in_ball_approx
    };

    cliques::find_all(cloud, in_ball, &mut cons);

    Filtration::new(cons.simplices)
}

/// Construct a Čech complex with exact smallest enclosing ball computation.
///
/// This is a convenience wrapper around `cech()` with `radius_tolerance = 0.0`.
/// Use this for guaranteed accuracy, though it may be slower for large cliques.
pub fn cech_exact<const N: usize>(
    space: &EuclideanCloud<N>,
    max_epsilon: Option<f64>,
    max_dim: Option<usize>,
) -> Filtration<Simplex>
{
    cech(space, max_epsilon, max_dim, 0.0)
}

fn in_ball_approx<const N: usize>(
    clique: &[usize],
    space: &EuclideanCloud<N>,
    cons: &Construction,
) -> Option<f64>
{
    let miniball = seb::larsson(clique, cons.tolerance, space);
    if miniball.r() > cons.max_epsilon {
        return None;
    }
    Some(miniball.r())
}

// TODO change algorithm used based on dimesion and size of clique
fn in_ball_exact<const N: usize>(
    clique: &[usize],
    space: &EuclideanCloud<N>,
    cons: &Construction,
) -> Option<f64>
{
    let clique_vec = clique.to_vec();
    let miniball = seb::welzl(&clique_vec, space);
    if miniball.r() > cons.max_epsilon {
        return None;
    }
    Some(miniball.r())
}
