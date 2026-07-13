


use crate::geometry::{ PointSet, seb };
use crate::topology::Filtration;
use crate::construction::{ cliques, Construction }; 


// NOTE without approximation should only be used on very sparse graphs as it gets
// extremely computationally expensive for larger cliques
pub fn cech(
    space : &PointSet,
    max_epsilon : Option<f64>,
    max_dim : Option<usize>,
    radius_tolerance : f64
) -> Filtration {
    let max_epsilon = max_epsilon.unwrap_or(f64::MAX / 2.0);
    let max_dim = max_dim.unwrap_or(usize::MAX - 1);
    let radius_tolerance = radius_tolerance.abs();

    let mut cons = Construction::new(max_dim, max_epsilon, radius_tolerance);
    cons.traverse_edges(space, 2.0 * max_epsilon, false);
    if cons.no_adjacency() { return Filtration::new(cons.simplices); }

    let in_ball = if radius_tolerance == 0.0 { in_ball_exact } else { in_ball_approx };

    let candidates : Vec<usize> = (0..space.len()).collect();
    cliques::bron_kerbosch(candidates, space, in_ball, &mut cons);

    cons.sort_simplices();
    Filtration::new(cons.simplices)
}


pub fn cech_exact(space : &PointSet, max_epsilon : Option<f64>, max_dim : Option<usize>) -> Filtration {
    cech(space, max_epsilon, max_dim, 0.0)
}


fn in_ball_approx(clique : &Vec<usize>, space : &PointSet, cons : &Construction) -> Option<f64> {
    let miniball = seb::larsson(clique, cons.tolerance, space);
    if miniball.r() > cons.max_epsilon { return None; }
    Some(miniball.r())
}


// TODO change algorithm used based on dimesion and size of clique
fn in_ball_exact(clique : &Vec<usize>, space : &PointSet, cons : &Construction) -> Option<f64> {
    let miniball = seb::welzl(clique, space);
    if miniball.r() > cons.max_epsilon { return None; }
    Some(miniball.r())
}
