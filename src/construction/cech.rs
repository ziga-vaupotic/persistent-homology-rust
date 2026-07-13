


use crate::geometry::{ PointSet, seb };
use crate::topology::Filtration;
use crate::construction::{ cliques, common::*, Distance };


// NOTE should only be used on very sparse graphs as it gets extremely computationally expensive for
// larger cliques
pub fn cech(
    space : &PointSet,
    max_epsilon : Option<f64>,
    max_dim : Option<usize>,
    radius_tolerance : f64 // radius <= (1 + radius_tolerance) r_optimal
) -> Filtration {
    let max_epsilon = max_epsilon.unwrap_or(f64::MAX / 2.0);
    let max_dim = max_dim.unwrap_or(usize::MAX - 1);
    let radius_tolerance = radius_tolerance.abs(); // did not feel liking returning Err

    let (mut simplices, adjacency, distance) = initialize_simplices(space, 2.0 * max_epsilon, false);
    if adjacency.is_empty() { return Filtration::new(simplices) }

    let f = if radius_tolerance == 0.0 { in_ball_exact } else { in_ball_approx };

    // TODO degeneracy ordering
    let candidates : Vec<usize> = (0..space.len()).collect(); // has to be ordered
    cliques::bron_kerbosch(
        Vec::new(),
        candidates,
        max_dim + 1,
        max_epsilon,
        radius_tolerance,
        space,
        &adjacency,
        &distance,
        f,
        &mut simplices
    );

    sort_simplices(&mut simplices);

    Filtration::new(simplices)
}


pub fn cech_exact(space : &PointSet, max_epsilon : Option<f64>, max_dim : Option<usize>) -> Filtration {
    cech(space, max_epsilon, max_dim, 0.0)
}


fn in_ball_approx(
    clique : &Vec<usize>,
    epsilon : f64,
    tolerance : f64,
    _distance : &Distance,
    space : &PointSet
) -> Option<f64> {
    let miniball = seb::larsson(clique, tolerance, space);
    if miniball.r() > epsilon { return None; }
    Some(miniball.r())
}


// TODO change algorithm used based on dimesion and size of clique
fn in_ball_exact(
    clique : &Vec<usize>,
    epsilon : f64,
    _tolerance : f64,
    _distance : &Distance,
    space : &PointSet
) -> Option<f64> {
    let miniball = seb::welzl(clique, space);
    if miniball.r() > epsilon { return None; }
    Some(miniball.r())
}
