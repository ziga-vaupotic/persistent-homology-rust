#[derive(Debug, Clone)]
/// A simplex in an abstract simplicial complex.
///
/// The simplex stores its vertex indices and a filtration value used to determine
/// the ordering in a filtration.
pub struct Simplex {
    pub vertices: Vec<usize>,
    pub filtration_value: f64, // this can be implemented as a trait instead later on!
}

use std::hash::{Hash, Hasher};

impl PartialEq for Simplex {
    fn eq(&self, other: &Self) -> bool {
        self.vertices == other.vertices
    }
}

impl Eq for Simplex {}

/* This will be used for indexing when computing homology */
impl Hash for Simplex {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.vertices.hash(state);
    }
}

impl Simplex {
    /// Create a new simplex with sorted vertex indices.
    pub fn new(vertices: Vec<usize>, filtration_value: f64) -> Self {
        let mut verts = vertices;
        verts.sort();
        Self {
            vertices: verts,
            filtration_value,
        }
    }

    /// Return the topological dimension of the simplex.
    pub fn dim(&self) -> usize {
        self.vertices.len() - 1
    }

    /// Compute the oriented boundary of the simplex.
    ///
    /// Returns pairs `(coefficient, face)` where coefficient is `1` or `-1`
    /// coming from the alternating boundary formula.
    pub fn boundary(&self) -> Vec<(i32, Simplex)> {
        let mut result = Vec::new();
        let n = self.vertices.len();
        for i in 0..n {
            let mut verts = self.vertices.clone();
            verts.remove(i);
            if !verts.is_empty() {
                let coeff = if i % 2 == 0 { 1 } else { -1 };
                let sub_simplex = Simplex::new(verts, 0.0); // filtration value doesn't matter for boundary
                result.push((coeff, sub_simplex));
            }
        }
        result
    }
}

/// Abstract simplical complex made of simplexes
pub struct SimplicialComplex {
    pub simplices: Vec<Simplex>,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_simplex_new_sorts_vertices() {
        let simplex = Simplex::new(vec![2, 0, 1], 1.0);
        assert_eq!(simplex.vertices, vec![0, 1, 2]);
        assert_eq!(simplex.filtration_value, 1.0);
    }

    #[test]
    fn test_simplex_dimension() {
        let simplex_0d = Simplex::new(vec![0], 0.0);
        assert_eq!(simplex_0d.dim(), 0);

        let simplex_1d = Simplex::new(vec![0, 1], 0.0);
        assert_eq!(simplex_1d.dim(), 1);

        let simplex_2d = Simplex::new(vec![0, 1, 2], 0.0);
        assert_eq!(simplex_2d.dim(), 2);
    }

    #[test]
    fn test_boundary_0_simplex() {
        let simplex = Simplex::new(vec![0], 0.0);
        let bound = simplex.boundary();
        assert_eq!(bound.len(), 0);
    }

    #[test]
    fn test_boundary_1_simplex() {
        let simplex = Simplex::new(vec![0, 1], 0.0);
        let bound = simplex.boundary();
        assert_eq!(bound.len(), 2);

        assert_eq!(bound[0], (1, Simplex::new(vec![1], 0.0)));
        assert_eq!(bound[1], (-1, Simplex::new(vec![0], 0.0)));
    }

    #[test]
    fn test_boundary_2_simplex() {
        let simplex = Simplex::new(vec![0, 1, 2], 0.0);
        let bound = simplex.boundary();
        assert_eq!(bound.len(), 3);

        assert_eq!(bound[0], (1, Simplex::new(vec![1, 2], 0.0)));
        assert_eq!(bound[1], (-1, Simplex::new(vec![0, 2], 0.0)));
        assert_eq!(bound[2], (1, Simplex::new(vec![0, 1], 0.0)));
    }
}
