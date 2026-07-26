//! Algebraic operations and matrix representations used by the persistent homology pipeline.
//!
//! This module handles discretisation of boundary operators into sparse matrices over $\mathbb{Z}_2$,
//! matrix reduction to compute persistence data, and extraction of persistence diagrams.

pub mod discretisation;
pub mod matrices;
pub mod persistence;
