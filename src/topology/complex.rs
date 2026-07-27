use crate::topology::Cell;

pub struct CellComplex<C: Cell> {
    pub cells: Vec<C>,
}
