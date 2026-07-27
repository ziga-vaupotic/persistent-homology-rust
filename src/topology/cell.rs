pub trait Cell: Clone + Eq {
    fn dim(&self) -> usize;
    fn boundary(&self) -> Vec<(i32, Self)>;
}

pub trait FilteredCell: Cell {
    fn filtration_value(&self) -> f64;
}
