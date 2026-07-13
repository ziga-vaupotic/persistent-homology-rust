


use crate::geometry::{ Ball, Point, PointSet };

use std::collections::HashMap;


// K. Fischer, B. Gärtner, and M. Kutz, Fast Smallest Enclosing Ball Computation in High Dimensions, 2003
// https://doi.org/10.1007/978-3-540-39658-1_57
// code tested in the article can be found at https://github.com/hbf/miniball, reposity owned by K. Fischer
pub fn fischer(_points : &Vec<usize>, _space : &PointSet, _d : HashMap<(usize, usize), f64>) -> Ball {
    Ball::new(Point::new(Vec::new()), 0.0)
}
