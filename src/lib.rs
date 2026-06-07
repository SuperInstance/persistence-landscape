//! # persistence-landscape
//!
//! **Persistence landscapes** for statistical topological data analysis.
//!
//! A persistence landscape is a sequence of piecewise-linear functions
//! `λ₁ ≥ λ₂ ≥ …` derived from a persistence diagram. They live in a Banach
//! space, enabling statistical analysis (averaging, hypothesis testing, etc.).
//!
//! # Example
//!
//! ```
//! use persistence_landscape::{PersistencePair, PersistenceLandscape, LandscapeNorm};
//!
//! let pairs = vec![PersistencePair::new(0.0, 1.0), PersistencePair::new(0.5, 2.0)];
//! let landscape = PersistenceLandscape::from_pairs(&pairs);
//! let v = landscape.lambda_k(1, 0.75);
//! let norm = LandscapeNorm::lp(&landscape, 2, 2.0);
//! ```

use std::cmp::Ordering;

// ---------------------------------------------------------------------------
// PersistencePair
// ---------------------------------------------------------------------------

/// A single persistence pair (birth, death).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PersistencePair {
    pub birth: f64,
    pub death: f64,
}

impl PersistencePair {
    /// Create a new persistence pair.
    pub fn new(birth: f64, death: f64) -> Self {
        Self { birth, death }
    }

    /// Persistence (death − birth).
    pub fn persistence(&self) -> f64 {
        self.death - self.birth
    }

    /// Midpoint of the pair.
    pub fn midpoint(&self) -> f64 {
        (self.birth + self.death) / 2.0
    }

    /// Height (half the persistence).
    pub fn height(&self) -> f64 {
        self.persistence() / 2.0
    }
}

impl PartialOrd for PersistencePair {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Eq for PersistencePair {}

impl Ord for PersistencePair {
    fn cmp(&self, other: &Self) -> Ordering {
        // Sort by decreasing persistence, then by birth
        match other.persistence().partial_cmp(&self.persistence()) {
            Some(Ordering::Equal) | None => self
                .birth
                .partial_cmp(&other.birth)
                .unwrap_or(Ordering::Equal),
            Some(o) => o,
        }
    }
}

// ---------------------------------------------------------------------------
// LandscapeFunction (internal)
// ---------------------------------------------------------------------------

/// A single piecewise-linear "tent" function λ_k derived from one persistence pair.
///
/// λ(x) = min(x − birth, death − x)⁺
///
/// which forms a tent peaking at the midpoint with height = persistence/2.
fn tent_value(pair: &PersistencePair, x: f64) -> f64 {
    let left = x - pair.birth;
    let right = pair.death - x;
    let m = left.min(right);
    if m > 0.0 { m } else { 0.0 }
}

// ---------------------------------------------------------------------------
// PersistenceLandscape
// ---------------------------------------------------------------------------

/// A persistence landscape: a sequence of piecewise-linear functions
/// `λ₁ ≥ λ₂ ≥ …` obtained by arranging the tent functions in decreasing order.
#[derive(Debug, Clone)]
pub struct PersistenceLandscape {
    /// Persistence pairs sorted by decreasing persistence.
    pairs: Vec<PersistencePair>,
}

impl PersistenceLandscape {
    /// Build a persistence landscape from a set of persistence pairs.
    pub fn from_pairs(pairs: &[PersistencePair]) -> Self {
        let mut sorted: Vec<PersistencePair> = pairs.to_vec();
        sorted.sort();
        Self { pairs: sorted }
    }

    /// Evaluate the k-th landscape function at point x (1-indexed).
    pub fn lambda_k(&self, k: usize, x: f64) -> f64 {
        if k == 0 || k > self.pairs.len() {
            return 0.0;
        }
        // The k-th landscape is the k-th largest tent value at x.
        let mut values: Vec<f64> = self.pairs.iter().map(|p| tent_value(p, x)).collect();
        values.sort_by(|a, b| b.partial_cmp(a).unwrap_or(Ordering::Equal));
        values[k - 1]
    }

    /// Number of landscape layers (same as number of pairs).
    pub fn num_layers(&self) -> usize {
        self.pairs.len()
    }

    /// Access the underlying pairs.
    pub fn pairs(&self) -> &[PersistencePair] {
        &self.pairs
    }

    /// Evaluate all landscape functions at x, returning λ₁(x), λ₂(x), …
    pub fn evaluate_all(&self, x: f64) -> Vec<f64> {
        let n = self.pairs.len();
        let mut values: Vec<f64> = self.pairs.iter().map(|p| tent_value(p, x)).collect();
        values.sort_by(|a, b| b.partial_cmp(a).unwrap_or(Ordering::Equal));
        values.resize(n, 0.0);
        values
    }
}

// ---------------------------------------------------------------------------
// LandscapeNorm
// ---------------------------------------------------------------------------

/// Compute norms of persistence landscapes.
pub struct LandscapeNorm;

impl LandscapeNorm {
    /// Compute the L^p norm of the k-th landscape function over a discretised domain.
    ///
    /// `num_samples` controls the resolution of the numerical integration.
    pub fn lp(landscape: &PersistenceLandscape, k: usize, p: f64) -> f64 {
        if landscape.pairs.is_empty() {
            return 0.0;
        }
        let lo = landscape
            .pairs
            .iter()
            .map(|pp| pp.birth)
            .fold(f64::INFINITY, f64::min);
        let hi = landscape
            .pairs
            .iter()
            .map(|pp| pp.death)
            .fold(f64::NEG_INFINITY, f64::max);

        let n = 200usize;
        let dx = (hi - lo) / n as f64;
        let mut sum = 0.0;
        for i in 0..=n {
            let x = lo + i as f64 * dx;
            let v = landscape.lambda_k(k, x);
            sum += v.powf(p) * dx;
        }
        sum.powf(1.0 / p)
    }

    /// Compute the supremum (L^∞) norm of the k-th landscape function.
    pub fn sup_norm(landscape: &PersistenceLandscape, k: usize) -> f64 {
        if k == 0 || k > landscape.pairs.len() {
            return 0.0;
        }
        // The k-th largest tent height
        landscape.pairs.get(k - 1).map(|p| p.height()).unwrap_or(0.0)
    }
}

// ---------------------------------------------------------------------------
// LandscapeAverage
// ---------------------------------------------------------------------------

/// Compute the average of multiple persistence landscapes.
pub struct LandscapeAverage;

impl LandscapeAverage {
    /// Average several landscapes by averaging their pair sets.
    ///
    /// For simplicity, this returns a landscape whose pairs are the
    /// pointwise-average birth/death of the corresponding pairs sorted by
    /// decreasing persistence. If the number of pairs differs, missing
    /// pairs are treated as (0,0) (zero persistence).
    pub fn average(landscapes: &[PersistenceLandscape]) -> PersistenceLandscape {
        if landscapes.is_empty() {
            return PersistenceLandscape { pairs: vec![] };
        }
        let max_len = landscapes.iter().map(|l| l.pairs.len()).max().unwrap_or(0);
        let n = landscapes.len() as f64;
        let mut avg_pairs = Vec::new();
        for i in 0..max_len {
            let mut b_sum = 0.0;
            let mut d_sum = 0.0;
            for l in landscapes {
                if let Some(p) = l.pairs.get(i) {
                    b_sum += p.birth;
                    d_sum += p.death;
                }
            }
            avg_pairs.push(PersistencePair::new(b_sum / n, d_sum / n));
        }
        PersistenceLandscape { pairs: avg_pairs }
    }
}

// ===========================================================================
// Tests
// ===========================================================================

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_pair_persistence() {
        let p = PersistencePair::new(1.0, 4.0);
        assert!((p.persistence() - 3.0).abs() < 1e-10);
        assert!((p.midpoint() - 2.5).abs() < 1e-10);
        assert!((p.height() - 1.5).abs() < 1e-10);
    }

    #[test]
    fn test_pair_ordering() {
        let p1 = PersistencePair::new(0.0, 5.0);
        let p2 = PersistencePair::new(0.0, 2.0);
        assert!(p1 < p2); // higher persistence comes first
    }

    #[test]
    fn test_tent_value_peak() {
        let p = PersistencePair::new(0.0, 2.0);
        assert!((tent_value(&p, 1.0) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_tent_value_zero_outside() {
        let p = PersistencePair::new(0.0, 2.0);
        assert!(tent_value(&p, -1.0).abs() < 1e-10);
        assert!(tent_value(&p, 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_lambda_k_first() {
        let pairs = vec![PersistencePair::new(0.0, 4.0)];
        let pl = PersistenceLandscape::from_pairs(&pairs);
        assert!((pl.lambda_k(1, 2.0) - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_lambda_k_beyond_returns_zero() {
        let pairs = vec![PersistencePair::new(0.0, 4.0)];
        let pl = PersistenceLandscape::from_pairs(&pairs);
        assert!((pl.lambda_k(2, 2.0)).abs() < 1e-10);
    }

    #[test]
    fn test_lambda_k_multiple_pairs() {
        let pairs = vec![
            PersistencePair::new(0.0, 4.0),
            PersistencePair::new(1.0, 3.0),
        ];
        let pl = PersistenceLandscape::from_pairs(&pairs);
        // λ₁ at x=2 should be max of tents = 2.0
        assert!((pl.lambda_k(1, 2.0) - 2.0).abs() < 1e-10);
        // λ₂ at x=2 should be second largest = 1.0
        assert!((pl.lambda_k(2, 2.0) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_num_layers() {
        let pairs = vec![
            PersistencePair::new(0.0, 1.0),
            PersistencePair::new(2.0, 5.0),
            PersistencePair::new(1.0, 3.0),
        ];
        let pl = PersistenceLandscape::from_pairs(&pairs);
        assert_eq!(pl.num_layers(), 3);
    }

    #[test]
    fn test_lp_norm_nonzero() {
        let pairs = vec![PersistencePair::new(0.0, 2.0)];
        let pl = PersistenceLandscape::from_pairs(&pairs);
        let norm = LandscapeNorm::lp(&pl, 1, 2.0);
        assert!(norm > 0.0);
    }

    #[test]
    fn test_sup_norm() {
        let pairs = vec![PersistencePair::new(0.0, 4.0)];
        let pl = PersistenceLandscape::from_pairs(&pairs);
        let sn = LandscapeNorm::sup_norm(&pl, 1);
        assert!((sn - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_sup_norm_beyond() {
        let pairs = vec![PersistencePair::new(0.0, 4.0)];
        let pl = PersistenceLandscape::from_pairs(&pairs);
        assert!((LandscapeNorm::sup_norm(&pl, 2)).abs() < 1e-10);
    }

    #[test]
    fn test_average_two_landscapes() {
        let l1 = PersistenceLandscape::from_pairs(&[PersistencePair::new(0.0, 2.0)]);
        let l2 = PersistenceLandscape::from_pairs(&[PersistencePair::new(2.0, 6.0)]);
        let avg = LandscapeAverage::average(&[l1, l2]);
        assert_eq!(avg.pairs.len(), 1);
        assert!((avg.pairs[0].birth - 1.0).abs() < 1e-10);
        assert!((avg.pairs[0].death - 4.0).abs() < 1e-10);
    }

    #[test]
    fn test_average_empty() {
        let avg = LandscapeAverage::average(&[]);
        assert_eq!(avg.num_layers(), 0);
    }
}
