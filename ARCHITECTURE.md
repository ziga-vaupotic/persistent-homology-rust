# Architecture Guide

This document describes the architecture and module structure of the persistent homology library.

## Overview

The library computes persistent homology of simplicial complexes, specifically focusing on Vietoris-Rips and Čech complexes constructed from point clouds. The computation pipeline follows these stages:

1. **Geometry**: Define a point cloud with an associated metric/inner product space
2. **Construction**: Build a filtered simplicial complex (Vietoris-Rips or Čech)
3. **Discretisation**: Convert the filtration into algebraic boundary matrices
4. **Algebra**: Reduce boundary matrices over $\mathbb{Z}_2$ and compute pivots
5. **Output**: Extract persistence diagrams and Betti numbers

## Module Structure

### `geometry` — Metric spaces and point clouds

- **`Point`**: A point in arbitrary-dimensional Euclidean space, backed by a dynamic vector.
- **`PointCloud<M>`**: A collection of points with an associated metric or inner product space.
- **Traits**:
  - `Metric`: Distance function satisfying metric axioms.
  - `Norm`: Norm (length) function; every norm induces a metric.
  - `InnerProduct`: Dot product; every inner product induces a norm.
  - `Euclidean`: Marker trait for Euclidean inner product spaces.
- **`EuclideanInnerProduct`**: Standard Euclidean dot product.
- **`Ball`**: Closed ball defined by center and radius; used for SEB tests in Čech complexes.
- **`seb` module**: Algorithms for computing smallest enclosing balls (Welzl, Larsson, Gärtner).

### `topology` — Simplicial complexes and filtrations

- **`Simplex`**: Represents a $k$-simplex as a sorted list of vertices and a filtration value.
  - `boundary()`: Computes oriented boundary as a list of $(k-1)$-faces.
  - `dim()`: Topological dimension (vertices.len() - 1).
- **`SimplicialComplex`**: Collection of simplices (not necessarily filtered).
- **`Filtration`**: An ordered sequence of simplices sorted by filtration value.
  - `complex_at(epsilon)`: Extract the simplicial complex at a given filtration level.
  - `simplices_of_dim(d)`: All simplices of dimension $d$.

### `construction` — Complex builders

- **`vietoris_rips<M>(space, max_epsilon, max_dim)`**: Construct a Vietoris-Rips complex.
  - Includes all $k$-simplices where all pairwise distances are $\leq \epsilon$.
  - Filtration value = maximum pairwise distance.
  - Works with any metric space.

- **`cech<M>(space, max_epsilon, max_dim, radius_tolerance)`**: Construct a Čech complex.
  - Includes a $k$-simplex if its smallest enclosing ball has radius $\leq \epsilon$.
  - Filtration value = SEB radius.
  - Requires Euclidean inner product space.
  - `radius_tolerance = 0.0` for exact computation (Welzl algorithm).
  - `radius_tolerance > 0` for approximate computation (Larsson algorithm, faster).

- **`cech_exact<M>(space, max_epsilon, max_dim)`**: Convenience wrapper for exact Čech.

Internal structures:
- **`Construction`**: State machine managing simplex generation and distance caching.
- **`cliques`**: Enumerate all $k$-cliques meeting geometric criteria.
- **`builder`**: Low-level construction utilities.

### `algebra` — Matrix operations and persistence

#### `matrices` — Boundary matrices and reduction

- **`BoundaryMatrix`**: Sparse column-form representation of a boundary operator over $\mathbb{Z}_2$.
  - Each column is a sorted list of row indices (entries are all 1 mod 2).
  - `column_indices`: Maps each column back to its global filtration index.
  - `reduce()`: Compute pivots using standard persistent homology algorithm.
  - `rank()`: Number of nonzero pivots.

- **`ReducedBoundaryMatrix`**: Result of reduction.
  - `matrix`: The reduced boundary matrix.
  - `low`: Pivot (low) indices; `low[j]` = largest row index in reduced column $j$, or `None` if zero.

- **Type aliases**:
  - `BoundaryMatrices = Vec<BoundaryMatrix>`
  - `ReducedBoundaryMatrices = Vec<ReducedBoundaryMatrix>`

#### `discretisation` — Filtration → matrices

- **`build_boundary_matrices(filtration)`**: Convert a filtration into boundary matrices.
  - One matrix per dimension $d \geq 1$.
  - Columns = boundaries of $d$-simplices (ordered by global filtration index).
  - Works over $\mathbb{Z}_2$ (all entries 0 or 1).

- **`reduce_boundary_matrices(matrices)`**: Reduce each matrix to persistent form.
  - Computes pivot (low) indices for extracting persistence pairs.
  - Uses Gaussian elimination with column operations over $\mathbb{Z}_2$.

#### `persistence` — Persistence pairs and diagrams

- **`PersistencePair`**: A single persistence interval.
  - `dimension`: Homology dimension.
  - `birth`: Global filtration index of birth simplex.
  - `death`: Global filtration index of death simplex (`None` = persists forever).
  - `is_finite()`, `is_infinite()`: Convenience methods.

- **`PersistenceDiagram`**: Collection of all persistence pairs.
  - `betti_at(epsilon)`: Compute Betti number (rank of homology) at filtration level.

- **Functions**:
  - `compute_persistence(matrix, dimension)`: Extract pairs from a single reduced matrix.
  - `compute_persistence_diagram(matrices)`: Combine pairs from all dimensions.

### `io` — File I/O

- **`import_point_cloud_csv<D, M>(path, geometry)`**: Load a point cloud from CSV.
  - Generic parameter `D`: expected coordinate dimension.
  - Validates consistent dimensions across all points.

- **`export_filtration_csv(path, filtration)`**: Write filtration to CSV.
  - Format: one simplex per row with filtration value, dimension, vertices.

## Typical Workflow

```rust
use persistent_homology::{
    geometry::{Point, PointCloud, EuclideanInnerProduct},
    construction::vietoris_rips,
    algebra::discretisation::{build_boundary_matrices, reduce_boundary_matrices},
    algebra::persistence::compute_persistence_diagram,
};

// 1. Load or create a point cloud
let points = vec![Point::new(vec![0.0]), Point::new(vec![1.0])];
let cloud = PointCloud::new(points, EuclideanInnerProduct)?;

// 2. Build a filtered complex
let filtration = vietoris_rips(&cloud, Some(2.0), Some(2));

// 3. Discretise into matrices
let matrices = build_boundary_matrices(&filtration);

// 4. Reduce matrices
let reduced = reduce_boundary_matrices(&matrices);

// 5. Extract persistence diagram
let diagram = compute_persistence_diagram(&reduced);

// 6. Analyze persistence
for pair in &diagram.pairs {
    println!("H_{}: [{}, {})", 
        pair.dimension, pair.birth, 
        pair.death.map_or("∞".to_string(), |d| d.to_string())
    );
}
```

## Key Design Decisions

1. **Sparse matrices over $\mathbb{Z}_2$**: Boundary matrices are stored in column-sparse form with all entries implicitly in $\{0, 1\}$ mod 2, reducing memory overhead.

2. **Filtration indices as coordinates**: Persistence pairs use global filtration indices (not filtration values) to avoid floating-point comparison issues and enable exact persistence computation.

3. **Trait-based geometry**: Point clouds are parameterized by a geometry type `M` (implementing `Metric`, `Norm`, or `InnerProduct`), enabling reuse across different metric spaces without runtime dispatch.

4. **Lazy construction**: Complex builders construct simplices on-demand using clique enumeration rather than storing all candidates, improving memory efficiency.

5. **$\mathbb{Z}_2$ reduction**: All matrix computations use $\mathbb{Z}_2$ coefficients (boolean algebra), simplifying the algorithm and improving performance.

## Performance Considerations

- **Clique enumeration**: Finding all maximal cliques is NP-hard; complexity grows exponentially with clique size. Limit `max_dim` for large point clouds.
- **SEB computation**: Welzl algorithm (exact Čech) is optimal but slower than Larsson (approximate). Use tolerance for speed.
- **Matrix reduction**: $O(n^3)$ in the worst case, where $n$ is the number of simplices. Sparse representation helps in practice.
- **Memory**: Filtrations can grow very large; consider streaming or sampling for massive datasets.
