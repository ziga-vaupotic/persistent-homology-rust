


use crate::geometry::{ Ball, Point, PointSet };


// B. Gärtner, Fast And Robust Smallest Enclosing Balls, 1999
// https://people.inf.ethz.ch/gaertner/subdir/software/miniball.html
// hits a bottle neck at higher dimensions as uses algorithm by Welzl as its base
// NOTE not fun to implement maybe at a later time
pub fn gartner(_P : &Vec<usize>, _space : &PointSet) -> Ball {
    Ball::new(Point::new(Vec::new()), 0.0)
}
