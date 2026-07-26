## Persistent homology in Rust
[![Rust](https://github.com/ziga-vaupotic/persistent-homology-rust/actions/workflows/rust.yml/badge.svg)](https://github.com/ziga-vaupotic/persistent-homology-rust/actions/workflows/rust.yml)
[![codecov](https://codecov.io/gh/ziga-vaupotic/persistent-homology-rust/graph/badge.svg?token=E3KUJJEBBC)](https://codecov.io/gh/ziga-vaupotic/persistent-homology-rust)
[![Crates.io](https://img.shields.io/crates/v/persistent-homology.svg)](https://crates.io/crates/persistent-homology)
[![Documentation](https://docs.rs/persistent-homology/badge.svg)](https://docs.rs/persistent-homology)
[![License](https://img.shields.io/crates/l/persistent-homology.svg)](LICENSE)

A lightweight Rust library for topological data analysis and persistent homology.

## Overview


`persistent-homology` serves as a lightweight implementation of persistent homology in Rust. The library supports abstractions such as simplicial complexes and allows for computing filtrations, and extracting topological information from point clouds over metric spaces.

The library is designed around generic abstractions for metric, normed, inner-product, and Euclidean spaces over $\mathbb R^n$. It supports a few basic algorithms for computing persistence.

This project was developed as part of the **Programiranje 2** course at the Faculty of Mathematics and Physics, University of Ljubljana.


## Features
The project currently supports topological data analysis over $(\mathbb R^d, d)$, where $d: \mathbb R^n \times \mathbb R^n \to \mathbb R$ is an arbitrary metric. The following core features are implemented:

- **Topology**: topological types and abstraction for persistence
    - Abstract simplex and simplicial complexes
    - Filtration over simplicial complexes
    - Abstraction for Metric, Normed, Inner-product and Euclidean spaces
- **Geometry**:
    - Points and point clouds in $(\mathbb R^n, d)$
    - Smallest enclosing balls for Euclidean spaces with several well-known algorithms
- **Filtration**: computation of filtration over point clouds
    - Clique based algorithms for persistence based on Bron-Kerbosch algorithm
    - Vietoris-Rips complex for arbitrary metric $d$
    - Čech complex (including approximations) for Euclidean space
- **Algebra**: matrix algebra needed for persistence:
    - Discretisation of boundary operators $\partial_{k}$ in simplicial complexes to matrices over the Galois field $GF(2)$
    - Reduction via optimised Gaussian elimination over $GF(2)$
- **Persistence**:
    - Computation of persistence pairs and persistence diagrams from reduced boundary matrices
    - Computation of Betti numbers

## Installation
To install the package simply run
```bash
cargo add persistent-homology
```
and add the following dependency to your `Cargo.toml`:
```toml
[dependencies] 
persistent-homology = "0.1.0"
```
For the latest version, see the package on [crates.io](https://crates.io/crates/persistent-homology).

## Example

The following example demonstrates the general workflow:

```rust
use persistent_homology::algebra::discretisation::{
    build_boundary_matrices, reduce_boundary_matrices,
};
use persistent_homology::algebra::persistence::compute_persistence_diagram;
use persistent_homology::construction::vietoris_rips;
use persistent_homology::geometry::EuclideanInnerProduct;


fn main() {
    // create point cloud over topology M
    let pointset = ...;

    // build filtration
    let filtration = vietoris_rips(&pointset, Some(4.0), Some(3));

    // build boundary matrices
    let matrices = build_boundary_matrices(&filtration);

    // reduce boundary matrices
    let reduced_matrices = reduce_boundary_matrices(&matrices);

    // compute persistence
    let persistence = compute_persistence_diagram(&reduced_matrices);
}
```

More complete examples can be found in the `examples` folder. For further documentation visit [docs.rs](https://docs.rs/persistent-homology/latest/persistent_homology/).

## Visual example
![image](logo.png)
This is a picture of Vietoris-Rips complex in Euclidean space $\mathbb{R}^2$ over the filtration with $\varepsilon=2.4$. The point cloud represents a noisy spiral.


## Use of Artificial Intelligence
Artificial intelligence was USED during the development of the library, however, the use was limited to writing (some) test and documentation. The
latter was checked by the authors.