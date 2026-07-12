


use crate::topology::{ Simplex, SimplicialComplex };


pub struct Filtration {
    pub simplices: Vec<Simplex>,
}



impl Filtration {

    pub fn new(simplices : Vec<Simplex>) -> Self {
        Self { simplices : simplices }
    }


    pub fn size(&self) -> usize {
        self.simplices.len()
    }

    pub fn max_dim(&self) -> usize {
        self.simplices.iter().max_by_key(|x| x.dim()).unwrap().dim()
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
