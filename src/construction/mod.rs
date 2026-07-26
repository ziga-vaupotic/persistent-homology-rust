//! Construction of simplicial complexes from point clouds.
//!
//! This module provides builders for Vietoris-Rips and Čech complexes, including utilities
//! for constructing geometric complexes from point clouds with customizable filtration parameters.
//! The construction process uses clique enumeration to build simplices efficiently.

mod cech;
mod rips;

pub use self::cech::*;
pub use self::rips::*;

mod builder;
mod cliques;

use self::builder::*;
