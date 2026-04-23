#[derive(Debug, Clone, PartialEq)]
pub struct Simplex {
    pub vertices: Vec<usize>,
    pub filtration_value: f64,
}

impl Simplex {
    pub fn dimension(&self) -> usize {
        self.vertices.len() - 1
    }
}

pub struct Filtration {
    pub simplices: Vec<Simplex>,
}