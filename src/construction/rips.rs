


use crate::construction::{ cliques, Construction }; 
use crate::geometry::PointSet;
use crate::topology::Filtration;

use itertools::Itertools;


/*
// for graph density
fn binom(a : usize, b : usize) -> usize {
    if b == 0 { return 1; }
    (a * binom(a - 1, b - 1)) / b
}
*/


pub fn vietoris_rips(space : &PointSet, max_epsilon : Option<f64>, max_dim : Option<usize>) -> Filtration {
    let max_epsilon = max_epsilon.unwrap_or(f64::MAX);
    let max_dim = max_dim.unwrap_or(usize::MAX - 1);

    let mut cons = Construction::new(max_dim, max_epsilon, 0.0);
    cons.traverse_edges(space, max_epsilon, true);
    if cons.no_adjacency() { return Filtration::new(cons.simplices); }

    /*
    let len_choose_2 = binom(space.len(), 2);
    let f : f64 = num_edges as f64 / len_choose_2 as f64;
    println!("graph density {} / {} = {:?}", num_edges, len_choose_2, f);
    */

    let candidates : Vec<usize> = (0..space.len()).collect();
    cliques::bron_kerbosch(candidates, space, rips_radius, &mut cons);

    cons.sort_simplices();
    Filtration::new(cons.simplices)
}


fn rips_radius(clique : &Vec<usize>, _space : &PointSet, cons : &Construction) -> Option<f64> {
    let mut max_d = 0.0;
    for v in (0..clique.len()).combinations(2) {
        let d = cons.distance[&(clique[v[0]], clique[v[1]])];
        max_d = if d >= max_d { d } else { max_d };
    }
    Some(max_d)
}


#[cfg(test)]
mod tests {

    use super::*;


    #[test]
    fn _test_vietoris_rips_three_points() {
        use std::f64::consts::PI;
        let points = vec![
            Point::new(vec![f64::cos(0.0), f64::sin(0.0)]),
            Point::new(vec![f64::cos(2.0 * PI / 3.0), f64::sin(2.0 * PI / 3.0)]),
            Point::new(vec![f64::cos(4.0 * PI / 3.0), f64::sin(4.0 * PI / 3.0)])
        ];
        let pointset = PointSet::new(points).expect("Pointset couldn't be generated.");

        let filtration = vietoris_rips(&pointset, None, None);

        // should have [0], [1], [2], [0, 1], [0, 2], [1, 2], [0, 1, 2]
        assert_eq!(filtration.simplices.len(), 7);

        // 0-simplices have filtration 0.0, 1 and 2-simplices 2 cos(PI / 6)
        for i in 0..3 {
            assert_eq!(filtration.simplices[i].filtration_value, 0.0);
        }
        for i in 3..7 {
            assert!((filtration.simplices[i].filtration_value - 2.0 * f64::cos(PI / 6.0)).abs() < 1e-12);
        }

        let dims: Vec<usize> = filtration.simplices.iter().map(|s| s.dim()).collect();
        assert_eq!(dims, vec![0, 0, 0, 1, 1, 1, 2]);

        assert_eq!(filtration.simplices[6].vertices, vec![0, 1, 2]);
    }

    #[test]
    fn test_vietoris_rips_single_point() {
        let points = vec![Point::new(vec![0.0, 0.0])];
        let pointset = PointSet::new(points).expect("Pointset couldn't be generated.");

        let filtration = vietoris_rips(&pointset, None, None);

        assert_eq!(filtration.simplices.len(), 1);
        assert_eq!(filtration.simplices[0].vertices, vec![0]);
        assert_eq!(filtration.simplices[0].filtration_value, 0.0);
    }
}
