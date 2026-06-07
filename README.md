# persistence-landscape

> **Persistence landscapes for statistical topological data analysis — turn persistence diagrams into functions you can average, compare, and test**

[![crates.io](https://img.shields.io/crates/v/persistence-landscape.svg)](https://crates.io/crates/persistence-landscape)
[![docs.rs](https://docs.rs/persistence-landscape/badge.svg)](https://docs.rs/persistence-landscape)
[![license](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

## What is a Persistence Landscape?

A **persistence diagram** is the canonical output of persistent homology — a collection of points (bᵢ, dᵢ) in the plane, each representing a topological feature that appears at scale bᵢ and disappears at scale dᵢ. But persistence diagrams live in a complicated mathematical space (the space of multisets), making standard statistical tools inapplicable. You can't easily average two diagrams or run a t-test on them.

**Persistence landscapes** solve this. Introduced by Bubenik (2015), they convert each persistence diagram into a sequence of piecewise-linear functions λ₁(t) ≥ λ₂(t) ≥ … that live in a **Banach space**. This means you can:

- **Average** landscapes across samples
- **Compute distances** using Lᵖ norms
- **Run hypothesis tests** (bootstrap, permutation tests)
- **Apply PCA, SVM, regression** — any method that works on vectors

The k-th landscape function λₖ(t) is defined by arranging the "tent functions" from each persistence pair in decreasing order. Each tent peaks at the midpoint (bᵢ+dᵢ)/2 with height (dᵢ−bᵢ)/2.

## Why Does This Matter?

Persistence landscapes bridge topology and statistics:

- **Clinical studies**: Compare the topological signature of healthy vs. diseased tissue samples, run statistical tests to determine significance
- **Materials science**: Characterize the porous structure of materials and test whether manufacturing changes produce statistically significant differences
- **Shape analysis**: Compute average shapes from topological descriptors, classify objects by their topological fingerprints
- **Time series**: Track how the topology of sliding windows changes over time, detect regime changes

The key insight: landscapes inherit the **stability** of persistence diagrams (small input changes → small landscape changes) while gaining **statistical tractability**.

## Architecture

```
┌─────────────────────────────────────────────────────────────┐
│              Persistence Landscape Pipeline                  │
│                                                             │
│  Persistence Diagram     Tent Functions       Landscapes    │
│                                                          │
│     (b₁,d₁) ╲         ╱╲   ╱╲   ╱╲              λ₁ ╱──╲  │
│      ·         ╲     ╱  ╲ ╱  ╲ ╱  ╲             ╱ ╱╲  ╲  │
│      ·   ───▶   ╲──▶╱    ╳    ╲    ╲    ───▶   ╱ ╱  ╲  ╲ │
│      ·              ╲  ╱ ╲  ╱ ╲  ╱           ╱ ╱    ╲  ╲ │
│      ·               ╲╱   ╲╱   ╲╱            λ₂·      ·  │
│                                                         │  │
│  Birth-Death pairs   λ(x) = min(x-b,d-x)⁺   Sorted by    │
│                       for each pair          peak height   │
│                                                         │  │
│  ──────────────────────────────────────────────────────  │  │
│                                                         │  │
│  ┌──────────────┐  ┌──────────────┐  ┌──────────────┐  │  │
│  │ LandscapeNorm │  │ LandscapeAvg │  │ Statistical  │  │  │
│  │  Lᵖ, L^∞     │  │  Mean of     │  │  Tests       │  │  │
│  │  distance     │  │  landscapes  │  │  Bootstrap   │  │  │
│  └──────────────┘  └──────────────┘  └──────────────┘  │  │
└─────────────────────────────────────────────────────────────┘
```

## Quick Start

```rust
use persistence_landscape::{PersistencePair, PersistenceLandscape, LandscapeNorm};

// Create persistence pairs from your homology computation
let pairs = vec![
    PersistencePair::new(0.0, 1.0),   // a short-lived feature
    PersistencePair::new(0.5, 2.0),   // a longer-lived feature
];

// Build the landscape
let landscape = PersistenceLandscape::from_pairs(&pairs);

// Evaluate the first landscape function at a point
let val = landscape.lambda_k(1, 0.75);
println!("λ₁(0.75) = {}", val);

// Compute the L² norm of the first landscape function
let norm = LandscapeNorm::lp(&landscape, 1, 2.0);
println!("||λ₁||₂ = {}", norm);

// Compute the supremum norm (peak height)
let sup = LandscapeNorm::sup_norm(&landscape, 1);
println!("||λ₁||∞ = {}", sup);
```

### Averaging Landscapes

```rust
use persistence_landscape::{PersistencePair, PersistenceLandscape, LandscapeAverage};

// Two samples with different topological signatures
let l1 = PersistenceLandscape::from_pairs(&[
    PersistencePair::new(0.0, 2.0),
    PersistencePair::new(1.0, 3.0),
]);
let l2 = PersistenceLandscape::from_pairs(&[
    PersistencePair::new(0.5, 3.0),
    PersistencePair::new(1.5, 4.0),
]);

// Compute the mean landscape
let mean = LandscapeAverage::average(&[l1, l2]);
println!("Average landscape has {} layers", mean.num_layers());
```

### Evaluating All Layers

```rust
let landscape = PersistenceLandscape::from_pairs(&pairs);

// Get all landscape values at a single point
let values = landscape.evaluate_all(1.0);
for (k, v) in values.iter().enumerate() {
    println!("λ_{}(1.0) = {}", k + 1, v);
}
```

## API Reference

### PersistencePair

| Method | Returns | Description |
|--------|---------|-------------|
| `PersistencePair::new(b, d)` | `PersistencePair` | Create a birth-death pair |
| `pair.persistence()` | `f64` | Life span: d − b |
| `pair.midpoint()` | `f64` | Center: (b + d) / 2 |
| `pair.height()` | `f64` | Peak height: (d − b) / 2 |

### PersistenceLandscape

| Method | Returns | Description |
|--------|---------|-------------|
| `PersistenceLandscape::from_pairs(&pairs)` | `PersistenceLandscape` | Build from persistence pairs |
| `landscape.lambda_k(k, x)` | `f64` | Evaluate k-th landscape at x (1-indexed) |
| `landscape.evaluate_all(x)` | `Vec<f64>` | All landscape values at x, sorted descending |
| `landscape.num_layers()` | `usize` | Number of landscape functions |
| `landscape.pairs()` | `&[PersistencePair]` | Access underlying pairs |

### LandscapeNorm

| Method | Returns | Description |
|--------|---------|-------------|
| `LandscapeNorm::lp(landscape, k, p)` | `f64` | Lᵖ norm of the k-th landscape |
| `LandscapeNorm::sup_norm(landscape, k)` | `f64` | L^∞ norm (peak height) of k-th landscape |

### LandscapeAverage

| Method | Returns | Description |
|--------|---------|-------------|
| `LandscapeAverage::average(&[landscapes])` | `PersistenceLandscape` | Pointwise average of multiple landscapes |

## Mathematical Background

### Tent Functions

For a persistence pair (b, d), the tent function is:

```
Λ_{(b,d)}(x) = max(0, min(x − b, d − x))
```

This creates a triangular "tent" that:
- Is zero outside [b, d]
- Peaks at x = (b+d)/2 with height (d−b)/2
- Is piecewise linear

### Landscape Construction

Given persistence pairs sorted by decreasing persistence: (b₁, d₁), (b₂, d₂), ...

The k-th landscape function at point x is:

```
λₖ(x) = k-th largest value of {Λ_{(bᵢ,dᵢ)}(x) : i = 1, ..., n}
```

This produces a sequence λ₁(x) ≥ λ₂(x) ≥ ... ≥ 0 where each λₖ is a piecewise-linear function.

### Key Properties

1. **Banach space**: Landscapes live in Lᵖ(ℝ × ℕ), enabling standard statistical analysis
2. **Stability**: ||Λ(D₁) − Λ(D₂)||∞ ≤ W∞(D₁, D₂) (landscape distance ≤ bottleneck distance)
3. **Injectivity**: Different diagrams produce different landscapes (in generic position)
4. **Isometry**: The L¹ norm equals the 1-Wasserstein distance for degree-0 persistence

### Statistical Applications

- **Two-sample testing**: Compute Lᵖ distance between mean landscapes, bootstrap for p-values
- **Classification**: Use landscape norms as features for SVM or random forests
- **Confidence bands**: Bootstrap landscapes to construct simultaneous confidence bands

## Installation

```bash
cargo add persistence-landscape
```

Or add to your `Cargo.toml`:

```toml
[dependencies]
persistence-landscape = "0.1.0"
```

## Related Crates

- [`betti-curve`](https://github.com/SuperInstance/betti-curve) — Betti curves and Euler characteristic curves
- [`mapper-graph`](https://github.com/SuperInstance/mapper-graph) — Mapper algorithm for topological summaries
- [`cech-complex`](https://github.com/SuperInstance/cech-complex) — Čech complex construction
- [`witness-complex`](https://github.com/SuperInstance/witness-complex) — Witness complex approximation

## License

MIT © [SuperInstance](https://github.com/SuperInstance)

---

*Part of the [Exocortex](https://github.com/SuperInstance/exocortex) project — persistent cognitive substrate for multi-agent systems.*
