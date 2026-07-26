#![allow(non_snake_case)]
//! A lightweight persistent homology library.
//!
//! This crate provides tools for building filtrations, constructing simplicial complexes,
//! discretising boundary operators into sparse matrices, and reducing boundary matrices
//! over the field `\mathbb{Z}_2` to compute persistent homology.

pub mod algebra;
pub mod construction;
pub mod geometry;
pub mod io;
pub mod topology;
