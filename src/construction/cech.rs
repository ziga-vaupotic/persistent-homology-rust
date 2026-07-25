


use crate::geometry::{ Euclidean, PointCloud, seb };
use crate::topology::Filtration;
use crate::construction::{ cliques, Construction }; 


pub fn cech<M>(
    space : &PointCloud<M>,
    max_epsilon : Option<f64>,
    max_dim : Option<usize>,
    radius_tolerance : f64
) -> Filtration
where
    M : Euclidean
{
    let max_epsilon = max_epsilon.unwrap_or(f64::MAX / 2.0);
    let max_dim = max_dim.unwrap_or(usize::MAX - 1);
    let radius_tolerance = radius_tolerance.abs();

    let mut cons = Construction::new(max_dim, max_epsilon, radius_tolerance);
    if !cons.traverse_edges(space, 2.0, false) { return Filtration::new(cons.simplices); }

    let in_ball = if radius_tolerance == 0.0 { in_ball_exact } else { in_ball_approx };

    cliques::find_all(space, in_ball, &mut cons);

    cons.sort_simplices();
    Filtration::new(cons.simplices)
}


pub fn cech_exact<M>(space : &PointCloud<M>, max_epsilon : Option<f64>, max_dim : Option<usize>) -> Filtration
where
    M : Euclidean
{
    cech(space, max_epsilon, max_dim, 0.0)
}


fn in_ball_approx<M>(clique : &Vec<usize>, space : &PointCloud<M>, cons : &Construction) -> Option<f64>
where
    M : Euclidean
{
    let miniball = seb::larsson(clique, cons.tolerance, space);
    if miniball.r() > cons.max_epsilon { return None; }
    Some(miniball.r())
}


// TODO change algorithm used based on dimesion and size of clique
fn in_ball_exact<M>(clique : &Vec<usize>, space : &PointCloud<M>, cons : &Construction) -> Option<f64>
where
    M : Euclidean
{
    let miniball = seb::welzl(clique, space);
    if miniball.r() > cons.max_epsilon { return None; }
    Some(miniball.r())
}
