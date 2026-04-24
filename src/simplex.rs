#[derive(Debug, Clone, PartialEq)]
pub struct Simplex {
    pub vertices: Vec<usize>,
    pub filtration_value: f64,
}

impl Simplex {
    pub fn new(vertices: Vec<usize>, filtration_value: f64) -> Self {
        let mut verts = vertices;
        verts.sort();
        Self {
            vertices: verts,
            filtration_value,
        }
    }

    pub fn dimension(&self) -> usize {
        self.vertices.len() - 1
    }
}

pub struct Filtration {
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
        assert_eq!(simplex_0d.dimension(), 0);

        let simplex_1d = Simplex::new(vec![0, 1], 0.0);
        assert_eq!(simplex_1d.dimension(), 1);

        let simplex_2d = Simplex::new(vec![0, 1, 2], 0.0);
        assert_eq!(simplex_2d.dimension(), 2);
    }
}