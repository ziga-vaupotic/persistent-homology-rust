


use crate::geometry::{ Ball, Point, PointSet };


// B. Gärtner, Fast And Robust Smallest Enclosing Balls, 1999
// https://people.inf.ethz.ch/gaertner/subdir/software/miniball.html
// NOTE hits a bottle neck at higher dimensions as uses welzl algorithm as its base
// NOTE not fun to implement as given in the paper
pub fn gartner(_P : &Vec<usize>, _space : &PointSet) -> Ball {
    Ball::new(Point::new(Vec::new()), 0.0)
}
