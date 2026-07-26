use nalgebra::DVector;

use crate::geometry::{Ball, Metric, Norm};

/// A point in an abstract Euclidean coordinate space.
///
/// Coordinates are stored in a dynamic vector so points can live in arbitrary dimension.
#[derive(Clone)]
pub struct Point {
    pub coords: DVector<f64>,
}

impl Point {
    /// Create a new point from a coordinate vector.
    pub fn new(coords: impl Into<DVector<f64>>) -> Self {
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
    pub fn standard_unit(i: usize, dim: usize) -> Self {
        let mut coords = vec![0.0; dim];
        coords[i] = 1.0;
        Point::new(coords)
    }
}

impl std::ops::Add<&Point> for &Point {
    type Output = Point;

    fn add(self, right: &Point) -> Point {
        let mut sum = Point::new(Vec::new());
        sum.coords = &self.coords + &right.coords;
        sum
    }
}

impl std::ops::Sub<&Point> for &Point {
    type Output = Point;

    fn sub(self, right: &Point) -> Point {
        let mut diff = Point::new(Vec::new());
        diff.coords = &self.coords - &right.coords;
        diff
    }
}

/// A collection of points with an associated metric or inner product space.
///
/// `PointCloud` enforces consistent point dimension and underlying geometry
/// via traits such as `Metric` and `InnerProduct`.
pub struct PointCloud<M> {
    // over R^n
    points: Vec<Point>,
    geometry: M, // Metric space e.g.
    dim: usize,  // This is added as an enforcement of consistency in a nodeset
}

impl<M> PointCloud<M>
where
    M: Copy,
{
    /// Create a point cloud with the given geometry.
    ///
    /// Returns an error when the input points do not all share the same coordinate dimension.
    pub fn new(points: Vec<Point>, geometry: M) -> Result<Self, String> {
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
    pub fn get(&self, i: usize) -> &Point {
        &self.points[i]
    }

    /// Return the metric or inner product instance used by this point cloud.
    pub fn get_geometry(&self) -> M {
        self.geometry
    }
}

impl<M> PointCloud<M>
where
    M: Metric,
{
    /// Compute the metric distance between two points.
    pub fn distance(&self, a: &Point, b: &Point) -> f64 {
        self.geometry.distance(a, b)
    }

    /// Check whether a point lies inside a closed ball.
    pub fn contained_in_ball(&self, ball: &Ball, point: &Point) -> bool {
        if ball.o().is_empty() {
            return false;
        }
        self.distance(ball.o(), point) <= ball.r()
    }
}

impl<M> PointCloud<M>
where
    M: Norm,
{
    /// Compute the squared norm of a point under the inner product.
    pub fn norm_squared(&self, a: &Point) -> f64 {
        self.geometry.norm(a) * self.geometry.norm(a)
    }

    /// Compute the Euclidean norm of a point under the inner product.
    pub fn norm(&self, a: &Point) -> f64 {
        self.geometry.norm(a)
    }
}
