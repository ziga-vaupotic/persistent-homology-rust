use nalgebra::SVector;

use crate::geometry::{Ball, Metric, Norm};
use std::ops::{Add, Sub};

/// A point in an abstract Euclidean coordinate space.
///
/// Coordinates are stored in a dynamic vector so points can live in arbitrary dimension.
#[derive(Clone)]
pub struct Point<const D: usize> {
    pub coords: SVector<f64, D>,
}

impl<const D: usize> Point<D> {
    /// Create a new point from a coordinate vector.
    pub fn new(coords: impl Into<SVector<f64, D>>) -> Self {
        Self {
            coords: coords.into(),
        }
    }

    /// Return the coordinate dimension of the point.
    pub fn len(&self) -> usize {
        self.coords.len()
    }

    /// Return `true` when the point has zero dimensions.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Return `true` when the point coordinate vector is numerically zero.
    pub fn is_zero(&self) -> bool {
        self.coords.iter().all(|x| x.abs() <= 1e-14)
    }

    /// Scale the point by a scalar multiplier.
    pub fn multiply(&mut self, lambda: f64) {
        self.coords = lambda * self.coords.clone();
    }

    /// Construct a standard basis vector in the given dimension.
    pub fn standard_unit(i: usize) -> Self {
        assert!(i < D);

        let mut coords = SVector::<f64, D>::zeros();
        coords[i] = 1.0;

        Point::new(coords)
    }
}

impl<const D: usize> Add<&Point<D>> for &Point<D> {
    type Output = Point<D>;

    #[inline]
    fn add(self, rhs: &Point<D>) -> Point<D> {
        Point {
            coords: self.coords + rhs.coords,
        }
    }
}

impl<const D: usize> Sub<&Point<D>> for &Point<D> {
    type Output = Point<D>;

    #[inline]
    fn sub(self, rhs: &Point<D>) -> Point<D> {
        Point {
            coords: self.coords - rhs.coords,
        }
    }
}

/// A collection of points with an associated metric or inner product space.
///
/// `PointCloud` enforces consistent point dimension and underlying geometry
/// via traits such as `Metric` and `InnerProduct`.
pub struct PointCloud<const D: usize, M> {
    // over R^n
    points: Vec<Point<D>>,
    geometry: M, // Metric space e.g.
    dim: usize,  // This is added as an enforcement of consistency in a nodeset
}

impl<const D: usize, M> PointCloud<D, M>
where
    M: Copy,
{
    /// Create a point cloud with the given geometry.
    ///
    /// Returns an error when the input points do not all share the same coordinate dimension.
    pub fn new(points: Vec<Point<D>>, geometry: M) -> Result<Self, String> {
        if points.is_empty() {
            return Ok(Self {
                points,
                geometry,
                dim: 0,
            });
        }

        let dim = points[0].coords.len();
        if !points.iter().all(|p| p.coords.len() == dim) {
            return Err("Inconsistent point dimensions".into());
        }

        Ok(Self {
            points,
            geometry,
            dim,
        })
    }

    /// Return the number of points in the cloud.
    pub fn len(&self) -> usize {
        self.points.len()
    }

    /// Return `true` when the point cloud contains no points.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Return the dimension of the underlying point space.
    pub fn dim(&self) -> usize {
        self.dim
    }

    /// Get a reference to the point at the given index.
    pub fn get(&self, i: usize) -> &Point<D> {
        &self.points[i]
    }

    /// Return the metric or inner product instance used by this point cloud.
    pub fn get_geometry(&self) -> M {
        self.geometry
    }
}

impl<const D: usize, M> PointCloud<D, M>
where
    M: Metric<D>,
{
    /// Compute the metric distance between two points.
    pub fn distance(&self, a: &Point<D>, b: &Point<D>) -> f64 {
        self.geometry.distance(a, b)
    }

    /// Check whether a point lies inside a closed ball.
    pub fn contained_in_ball(&self, ball: &Ball<D>, point: &Point<D>) -> bool {
        if ball.o().is_empty() {
            return false;
        }
        self.distance(ball.o(), point) <= ball.r()
    }
}

impl<const D: usize, M> PointCloud<D, M>
where
    M: Norm<D>,
{
    /// Compute the squared norm of a point under the inner product.
    pub fn norm_squared(&self, a: &Point<D>) -> f64 {
        self.geometry.norm(a) * self.geometry.norm(a)
    }

    /// Compute the Euclidean norm of a point under the inner product.
    pub fn norm(&self, a: &Point<D>) -> f64 {
        self.geometry.norm(a)
    }
}
