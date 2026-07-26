use crate::topology::{Simplex, SimplicialComplex};

/// A filtration is an ordered sequence of simplical complexes.
///
/// It stores every simplex separately for performance,
/// however, a simplical complex can be extracted.
pub struct Filtration {
    pub simplices: Vec<Simplex>,
}

impl Filtration {
    /// Create a new filtration from an ordered list of simplices.
    ///
    /// # Arguments
    ///
    /// * `simplices` - A list of (order) simplices.
    pub fn new(simplices: Vec<Simplex>) -> Self {
        Self { simplices }
    }

    /// Return the number of simplices in the filtration.
    pub fn len(&self) -> usize {
        self.simplices.len()
    }

    /// Return true if filtration is empty.
    pub fn is_empty(&self) -> bool {
        self.simplices.len() == 0
    }

    /// Return the maximal dimension of simplical complex.
    pub fn max_dim(&self) -> usize {
        self.simplices.iter().max_by_key(|x| x.dim()).unwrap().dim()
    }

    /// Return the filtration values of largest simplical complex.
    pub fn max_filtration_value(&self) -> f64 {
        self.simplices
            .iter()
            .max_by(|x, y| x.filtration_value.partial_cmp(&y.filtration_value).unwrap())
            .unwrap()
            .filtration_value
    }

    /// Return the simplices of dim of largest simplical complex.
    /// 
    /// # Arguments
    ///
    /// * `dim` - dimension of simplical complex.
    pub fn simplices_of_dim(&self, dim: usize) -> Vec<Simplex> {
        let mut simplices: Vec<Simplex> = Vec::new();
        for x in self.simplices.iter() {
            if x.dim() != dim {
                continue;
            }
            simplices.push(x.clone());
        }
        simplices
    }

    /// Return the simplicial complex consisting of all simplices with filtration
    /// value less than or equal to `epsilon`.
    /// 
    /// # Arguments
    ///
    /// * `epsilon` - 'epsilon' of simplical complex.
    pub fn complex_at(&self, epsilon: f64) -> SimplicialComplex {
        SimplicialComplex {
            simplices: self
                .simplices
                .iter()
                .filter(|s| s.filtration_value <= epsilon)
                .cloned()
                .collect(),
        }
    }
}
