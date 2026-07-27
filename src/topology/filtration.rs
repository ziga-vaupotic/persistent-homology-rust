use crate::topology::{CellComplex, FilteredCell};

/// A filtration is an ordered sequence of simplical complexes.
///
/// It stores every simplex separately for performance,
/// however, a simplical complex can be extracted.
pub struct Filtration<C: FilteredCell> {
    pub cells: Vec<C>,
}

impl<C: FilteredCell> Filtration<C> {
    /// Create a new filtration from an ordered list of simplices.
    ///
    /// # Arguments
    ///
    /// * `simplices` - A list of (order) simplices.
    pub fn new(mut cells: Vec<C>) -> Self {
        Self::sort_cells(&mut cells);
        Self { cells }
    }

    /// Return the number of simplices in the filtration.
    pub fn len(&self) -> usize {
        self.cells.len()
    }

    /// Return true if filtration is empty.
    pub fn is_empty(&self) -> bool {
        self.cells.len() == 0
    }

    /// Return the maximal dimension of simplical complex.
    pub fn max_dim(&self) -> usize {
        self.cells.iter().max_by_key(|x| x.dim()).unwrap().dim()
    }

    /// Return the filtration values of largest simplical complex.
    pub fn max_filtration_value(&self) -> f64 {
        self.cells
            .iter()
            .max_by(|x, y| {
                x.filtration_value()
                    .partial_cmp(&y.filtration_value())
                    .unwrap()
            })
            .unwrap()
            .filtration_value()
    }

    // Sort filtration, without return value.
    pub fn sort(&mut self) {
        Self::sort_cells(&mut self.cells);
    }

    // Sort a list of celles based on filtration value
    ///
    /// # Arguments
    ///
    /// * `cells` - list of a mutable cells.
    fn sort_cells(cells: &mut [C]) {
        cells.sort_by(|a, b| {
            a.filtration_value()
                .partial_cmp(&b.filtration_value())
                .unwrap()
                .then(a.dim().cmp(&b.dim()))
                .then(a.boundary().len().cmp(&b.boundary().len()))
        });
    }

    /// Return the simplices of dim of largest cells complex.
    ///
    /// # Arguments
    ///
    /// * `dim` - dimension of cells complex.
    pub fn cells_of_dim(&self, dim: usize) -> Vec<C> {
        let mut cells: Vec<C> = Vec::new();
        for x in self.cells.iter() {
            if x.dim() != dim {
                continue;
            }
            cells.push(x.clone());
        }
        cells
    }

    /// Return the simplicial complex consisting of all simplices with filtration
    /// value less than or equal to `epsilon`.
    ///
    /// # Arguments
    ///
    /// * `epsilon` - 'epsilon' of simplical complex.
    pub fn complex_at(&self, epsilon: f64) -> CellComplex<C> {
        CellComplex {
            cells: self
                .cells
                .iter()
                .filter(|s| (*s).filtration_value() <= epsilon)
                .cloned()
                .collect(),
        }
    }
}
