//! Geometry primitives and metric-space abstractions used by the filtration builders.
//!
//! This module provides point cloud representations, metric and inner product traits,
//! Euclidean vector operations, and balls for containment tests.

pub mod seb; // smallest enclosing ball

mod ball;
mod euclidean;
mod point_cloud;
mod space;

pub use self::ball::*;
pub use self::euclidean::*;
pub use self::point_cloud::*;
pub use self::space::*;
