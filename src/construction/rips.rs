


use crate::construction::{ cliques, Construction }; 
use crate::geometry::{ Metric, PointCloud };
use crate::topology::Filtration;

use itertools::Itertools;


pub fn vietoris_rips<M>(
    space : &PointCloud<M>,
    max_epsilon : Option<f64>,
    max_dim : Option<usize>
) -> Filtration
where
    M : Metric
{
    let max_epsilon = max_epsilon.unwrap_or(f64::MAX);
    let max_dim = max_dim.unwrap_or(usize::MAX - 1);

    let mut cons = Construction::new(max_dim, max_epsilon, 0.0);
    if !cons.traverse_edges(space, 1.0, true) { return Filtration::new(cons.simplices); }

    cliques::find_all(space, rips_radius, &mut cons);

    cons.sort_simplices();
    Filtration::new(cons.simplices)
}


fn rips_radius<M>(clique : &Vec<usize>, _space : &PointCloud<M>, cons : &Construction) -> Option<f64>
where
    M : Metric
{
    let mut max_d = 0.0;
    for v in (0..clique.len()).combinations(2) {
        let d = cons.distance[&(clique[v[0]], clique[v[1]])];
        max_d = if d >= max_d { d } else { max_d };
    }
    Some(max_d)
}
