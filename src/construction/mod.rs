


mod common;
mod cliques;

mod cech;
mod rips;

pub use self::cech::*;
pub use self::rips::*;

use std::collections::HashMap;
type Distance = HashMap<(usize, usize), f64>;
type Adjacency = HashMap<usize, Vec<usize>>;
