## Persistent homology in Rust
[![Rust](https://github.com/ziga-vaupotic/persistent-homology-rust/actions/workflows/rust.yml/badge.svg)](https://github.com/ziga-vaupotic/persistent-homology-rust/actions/workflows/rust.yml)
[![codecov](https://codecov.io/gh/ziga-vaupotic/persistent-homology-rust/graph/badge.svg?token=E3KUJJEBBC)](https://codecov.io/gh/ziga-vaupotic/persistent-homology-rust)
[![Crates.io](https://img.shields.io/crates/v/persistent-homology.svg)](https://crates.io/crates/persistent-homology)
[![Documentation](https://docs.rs/persistent-homology/badge.svg)](https://docs.rs/persistent-homology)
[![License](https://img.shields.io/crates/l/persistent-homology.svg)](LICENSE)

## Overview

This project servers as a lightweight implementation of persistent homology in Rust. The project was created as a part of the Programiranje 2 course at Faculty of Mathematics and Physics in Ljubljana.

Examples can be found in the `examples` folder.

![image](logo.png)
This is a picture of Vietoris-Rips complex in Euclidean space $\mathbb{R}^2$ over the filtration with $\varepsilon=2.4$. The point cloud represents a noisy spiral.
***

## Features
The project currently supports:
- Point clouds over arbitrary metric spaces
- Simplexes and simplicial complexes
- Filtration of abstract types
- clique based Vietoris-Rips and Čech complexes
- Boundary operator discretisation
- Boundary matrices and $\mathbb{Z}_2$ reduction
- Betti number computation


## On the use of AI
Artificial intelligence was USED during the development of the library, however, the use was limited to writing (some) test and documentation. The
latter was checked by the authors.