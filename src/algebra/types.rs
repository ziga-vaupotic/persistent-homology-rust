pub struct BoundaryMatrix {
    columns: Vec<Vec<usize>>,
}

impl BoundaryMatrix {
    pub fn new(columns: Vec<Vec<usize>>) -> Self {
        Self { columns }
    }

    pub fn columns(&self) -> &Vec<Vec<usize>> {
        &self.columns
    }
}
pub type BoundaryMatrices = Vec<BoundaryMatrix>;