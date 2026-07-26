//! Topological primitives for simplicial complexes and filtrations.
//!
//! This module exposes simplex definitions, simplicial complexes, and filtration
//! utilities used by the persistent homology pipeline.

mod filtration;
mod simplex;

pub use self::filtration::*;
pub use self::simplex::*;
