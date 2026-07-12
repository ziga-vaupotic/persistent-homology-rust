
use crate::algebra::PersistenceDiagram;

pub fn betti_at(
    diagram: &PersistenceDiagram,
    epsilon: usize
) -> usize {
    diagram.pairs
        .iter()
        .filter(|pair| {
            pair.birth <= epsilon &&
            pair.death.map_or(true, |death| epsilon < death)
        })
        .count()
}
