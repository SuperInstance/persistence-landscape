# persistence-landscape

> **Persistence landscapes for statistical topological data analysis**

[![crates.io](https://img.shields.io/crates/v/persistence-landscape.svg)](https://crates.io/crates/persistence-landscape)
[![license](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

Persistence landscapes transform persistence diagrams into functions, enabling statistical analysis. While persistence diagrams live in a non-vector space (making statistics hard), landscapes live in a Banach space — you can average them, compute norms, and run standard statistical tests.

## Key Idea

Given persistence pairs {(bᵢ, dᵢ)}, the landscape Λ(t) is a piecewise-linear function that peaks at each (bᵢ+dᵢ)/2 and goes to zero at bᵢ and dᵢ. Multiple landscape functions λₖ are ordered by peak height.

## Why Landscapes?

- **Vector space**: Can compute means, variances
- **Stability**: Small changes in diagram → small changes in landscape
- **Statistical tests**: Bootstrap, hypothesis testing on landscapes
- **Lp norms**: Natural distance metric between landscapes

## Installation

```toml
[dependencies]
persistence-landscape = "0.1.0"
```

## License

MIT © [SuperInstance](https://github.com/SuperInstance)

---

*Part of the [Exocortex](https://github.com/SuperInstance/exocortex) project — persistent cognitive substrate for multi-agent systems.*
