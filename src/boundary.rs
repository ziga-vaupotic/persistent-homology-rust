use crate::simplex::Simplex;

pub fn boundary(simplex: &Simplex) -> Vec<(i32, Simplex)> {
    let mut result = Vec::new();
    let n = simplex.vertices.len();
    for i in 0..n {
        let mut verts = simplex.vertices.clone();
        verts.remove(i);
        if !verts.is_empty() {
            let coeff = if i % 2 == 0 { 1 } else { -1 };
            let sub_simplex = Simplex::new(verts, 0.0); // filtration value doesn't matter for boundary
            result.push((coeff, sub_simplex));
        }
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_boundary_0_simplex() {
        let simplex = Simplex::new(vec![0], 0.0);
        let bound = boundary(&simplex);
        assert_eq!(bound.len(), 0);
    }

    #[test]
    fn test_boundary_1_simplex() {
        let simplex = Simplex::new(vec![0, 1], 0.0);
        let bound = boundary(&simplex);
        assert_eq!(bound.len(), 2);

        assert_eq!(bound[0], (1, Simplex::new(vec![1], 0.0)));
        assert_eq!(bound[1], (-1, Simplex::new(vec![0], 0.0)));
    }

    #[test]
    fn test_boundary_2_simplex() {
        let simplex = Simplex::new(vec![0, 1, 2], 0.0);
        let bound = boundary(&simplex);
        assert_eq!(bound.len(), 3);

        assert_eq!(bound[0], (1, Simplex::new(vec![1, 2], 0.0)));
        assert_eq!(bound[1], (-1, Simplex::new(vec![0, 2], 0.0)));
        assert_eq!(bound[2], (1, Simplex::new(vec![0, 1], 0.0)));
    }
}