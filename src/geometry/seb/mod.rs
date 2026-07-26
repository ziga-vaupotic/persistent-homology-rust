//! Smallest enclosing ball algorithms.
//!
//! This module provides multiple algorithms for computing the smallest enclosing ball (SEB)
//! of a set of points. The list of algorithms that are currently supported follows:
//! - `welzl`: Exact algorithm using recursive randomized approach (optimal but slower).
//! - `larsson`: Approximate algorithm using heuristics (faster, suitable for approximations).
//! - `fischer`: Alternative SEB algorithm variant.
//! - `gartner`: Alternative SEB algorithm variant.

mod fischer;
mod gartner;
mod larsson;
mod welzl;

pub use self::fischer::*;
pub use self::gartner::*;
pub use self::larsson::*;
pub use self::welzl::*;
