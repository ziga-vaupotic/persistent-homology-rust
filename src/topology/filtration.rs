


use crate::topology::simplex::Simplex;
use crate::topology::simplical_complex::SimplicialComplex;


pub struct Filtration {
    pub simplices: Vec<Simplex>,
}



impl Filtration {

    pub fn new(simplices : Vec<Simplex>) -> Self {
        Self { simplices : simplices }
    }


    // Instead of copying elements we could also get only a vector of references to simplices!
    pub fn get_simplicial_complex(&self, epsilon: f64) -> SimplicialComplex {
        SimplicialComplex {
            simplices: self.simplices
                .iter()
                .filter(|s| s.filtration_value <= epsilon)
                .cloned()
                .collect(),
        }
    }

}
