//! Persistence diagram computation for reduced boundary matrices.
//!
//! This module converts reduced boundary matrices into persistence diagram

use crate::algebra::matrices::{ReducedBoundaryMatrices, ReducedBoundaryMatrix};
use std::collections::HashSet;

/// A single persistence interval in a given dimension.
///
/// `birth` and `death` are global filtration indices of the corresponding simplices.
/// A `death` of `None` indicates the interval persists forever (born but never dies).
pub struct PersistencePair {
    pub dimension: usize,
    pub birth: usize,
    pub death: Option<usize>,
}

impl PersistencePair {
    /// Check if this pair represents a finite interval (birth and death).
    pub fn is_finite(&self) -> bool {
        self.death.is_some()
    }

    /// Check if this pair represents an infinite interval (born but never dies).
    pub fn is_infinite(&self) -> bool {
        self.death.is_none()
    }
}

/// A collection of persistence intervals.
///
/// The persistence diagram is stored as a vector of pairs, one for each finite or
/// infinite homology class.
pub struct PersistenceDiagram {
    pub pairs: Vec<PersistencePair>,
}

impl PersistenceDiagram {
    /// Compute the Betti number at the given filtration index.
    ///
    /// This counts intervals that are born on or before `epsilon` and not yet died.
    pub fn betti_at(&self, epsilon: usize) -> usize {
        self.pairs
            .iter()
            .filter(|pair| pair.birth <= epsilon && pair.death.is_none_or(|death| epsilon < death))
            .count()
    }
}

/// Compute persistence pairs for a single reduced boundary matrix.
///
/// The returned pairs represent finite intervals for each pivot (paired birth-death pairs)
/// and infinite intervals for zero columns whose births were never paired as deaths.
///
/// # Arguments
///
/// * `matrix` - A reduced boundary matrix with computed pivot indices.
/// * `dimension` - The homology dimension of this matrix.
///
/// # Returns
///
/// A vector of `PersistencePair` structs representing all persistence intervals at the given dimension.
pub fn compute_persistence(
    matrix: &ReducedBoundaryMatrix,
    dimension: usize,
) -> Vec<PersistencePair> {
    let mut pairs = Vec::new();

    // Birth simplices that get killed
    let mut paired_births = HashSet::new();

    // Finite intervals
    for (local_col, pivot_row) in matrix.low.iter().enumerate() {
        if let Some(birth) = pivot_row {
            let death = matrix.matrix.column_indices[local_col];

            paired_births.insert(*birth);

            pairs.push(PersistencePair {
                dimension,
                birth: *birth,
                death: Some(death),
            });
        }
    }

    // Infinite intervals:
    // zero reduced columns whose birth simplex
    // never appears as a pivot row.
    for (local_col, col) in matrix.matrix.columns.iter().enumerate() {
        if col.is_empty() {
            let birth = matrix.matrix.column_indices[local_col];

            if !paired_births.contains(&birth) {
                pairs.push(PersistencePair {
                    dimension,
                    birth,
                    death: None,
                });
            }
        }
    }

    pairs
}

/// Compute a complete persistence diagram from reduced boundary matrices across all dimensions.
///
/// Combines persistence pairs from each dimension into a single `PersistenceDiagram`.
/// The dimension index automatically corresponds to the homological dimension.
///
/// # Arguments
///
/// * `matrices` - A slice of reduced boundary matrices, one per homology dimension.
///
/// # Returns
///
/// A `PersistenceDiagram` containing all persistence intervals from all dimensions.
pub fn compute_persistence_diagram(matrices: &ReducedBoundaryMatrices) -> PersistenceDiagram {
    let mut pairs = Vec::new();

    for (dimension, matrix) in matrices.iter().enumerate() {
        pairs.extend(compute_persistence(matrix, dimension));
    }

    PersistenceDiagram { pairs }
}
