//! Geometry primitives and metric-space abstractions used by the filtration builders.
//!
//! This module provides point cloud representations, metric and inner product traits,
//! Euclidean vector operations, and balls for containment tests.

pub mod seb; // smallest enclosing ball

mod ball;
mod points;
mod space;

pub use self::ball::*;
pub use self::points::*;
pub use self::space::*;
