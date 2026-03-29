// Copyright (c) 2025 Shenghao Yang. All rights reserved.
// Licensed under the MIT License. See LICENSE-MIT for details.

use super::degree_dist::*;
use crate::types::ProbValue;

/// Sample a degree from a CDF.
/// r must be a uniformly distributed random number in [0, m), where 
/// m is the last value of the CDF.
pub fn sample_degree_from_cdf(cdf: &[ProbValue], r: usize) -> usize {
    let mut i = 0;
    while r >= cdf[i] as usize {
        i += 1;
    }
    i
}
/// Sample a degree from a PDF.
/// r must be a uniformly distributed random number in [0, m), where 
/// m is the sum of the PDF.
pub fn sample_degree_from_pdf(pdf: &[ProbValue], r: usize) -> usize {
    let cdf = pdf_to_cdf(pdf);
    sample_degree_from_cdf(&cdf, r)
}

/// Deterministic degree counts for the first \(n\) coded vectors (doc-scheme: “equal distance”
/// degree distribution sampling).
///
/// [`DeterministicDegreeValues::generate`] returns a multiset of \(n\) degrees (sorted ascending);
/// any permutation is valid for assignment to coded vector indices.
pub struct DeterministicDegreeValues {
    pdf: Vec<u32>,
    m: u32,
    n_d: Vec<usize>,
    last_n: u32,
}

impl DeterministicDegreeValues {
    /// Builds from a CDF in the same convention as the rest of this module (`cdf[0] == 0`,
    /// `cdf[dmax] == m`). `m` is taken as `cdf.last()`.
    pub fn new_from_cdf(cdf: Vec<u32>) -> Self {
        assert!(cdf.len() >= 2, "CDF must cover at least degree 1 (len >= 2)");
        let m = *cdf.last().expect("cdf non-empty");
        if !validate_cdf(&cdf, m as usize) {
            panic!("Invalid CDF: {:?}", cdf);
        }
        let pdf = cdf_to_pdf(&cdf);
        let len = pdf.len();
        Self {
            pdf,
            m,
            n_d: vec![0; len],
            last_n: 0,
        }
    }

    /// Same apportionment from an explicit PDF; `m` must equal `pdf.iter().sum()`.
    pub fn new_from_pdf(pdf: Vec<u32>, m: u32) -> Self {
        let sum: u64 = pdf.iter().map(|&x| x as u64).sum();
        assert_eq!(
            sum,
            m as u64,
            "PDF sum must equal m"
        );
        let len = pdf.len();
        Self {
            pdf,
            m,
            n_d: vec![0; len],
            last_n: 0,
        }
    }

    /// Degree counts `n_d` for the last `n` used in [`Self::recompute_counts`] / [`Self::generate`].
    #[inline]
    pub fn counts(&self) -> &[usize] {
        &self.n_d
    }


    /// Recomputes `n_d` for exactly `n` coded vectors (largest-remainder apportionment).
    pub fn recompute_counts(&mut self, n: u32) {
        let n_usize = n as usize;

        let mut counts = vec![0usize; self.pdf.len()];
        if n == 0 {
            self.n_d = counts;
            self.last_n = 0;
            return;
        }

        let mut sum_q = 0usize;
        // (remainder numer mod m, degree d) — tie-break smaller `d` first after sorting by rem desc
        let mut meta: Vec<(u32, usize)> = Vec::new();

        // With `n: u32` and `pdf[d]: u32`, the product `n * pdf[d]` always fits in `u64`
        // (max = (2^32-1)^2 < 2^64), so we can stay entirely in integer arithmetic.
        let m_u64 = self.m as u64;
        let n_u64 = n as u64;

        for d in 1..self.pdf.len() {
            let pd_u64 = self.pdf[d] as u64;

            let prod_u64 = n_u64 * pd_u64;
            let q = (prod_u64 / m_u64) as usize;
            let r = (prod_u64 % m_u64) as u32;
            counts[d] = q;
            sum_q = sum_q.saturating_add(q);
            meta.push((r, d));
        }

        meta.sort_by(|a, b| b.0.cmp(&a.0).then_with(|| a.1.cmp(&b.1)));

        let mut deficit = n_usize.saturating_sub(sum_q);
        debug_assert!(deficit <= meta.len() || meta.is_empty());
        let mut i = 0;
        while deficit > 0 && i < meta.len() {
            let d = meta[i].1;
            counts[d] += 1;
            deficit -= 1;
            i += 1;
        }
        assert_eq!(
            deficit, 0,
            "apportionment deficit should be zero for a valid PDF (sum pdf == m)"
        );

        self.n_d = counts;
        self.last_n = n;
    }

    /// Returns the multiset of `n` degree values (non-decreasing order: all degree-1, then degree-2, …).
    pub fn generate(&mut self, n: u32) -> Vec<usize> {
        self.recompute_counts(n);
        let mut seq = Vec::with_capacity(n as usize);
        for d in 1..self.n_d.len() {
            for _ in 0..self.n_d[d] {
                seq.push(d);
            }
        }
        debug_assert_eq!(seq.len(), n as usize);
        seq
    }

    /// For extending from `n` to `n_prime` coded vectors, returns the multiset of the additional
    /// `n_prime - n` degrees (`n'_d - n_d` expanded as values), in non-decreasing order.
    pub fn additional_degrees(&mut self, n: u32, n_prime: u32) -> Vec<usize> {
        assert!(
            n_prime >= n,
            "n_prime must be >= n"
        );
        self.recompute_counts(n);
        let old = self.n_d.clone();
        self.recompute_counts(n_prime);
        let mut seq = Vec::with_capacity((n_prime - n) as usize);
        for d in 1..self.n_d.len() {
            let add = self.n_d[d].saturating_sub(old[d]);
            for _ in 0..add {
                seq.push(d);
            }
        }
        debug_assert_eq!(seq.len(), (n_prime - n) as usize);
        seq
    }
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_deterministic_degree_values_uniform() {
        let pdf = vec![0u32, 1, 1, 1];
        let mut g = DeterministicDegreeValues::new_from_pdf(pdf, 3);
        let seq = g.generate(3);
        assert_eq!(seq, vec![1, 2, 3]);
    }

    #[test]
    fn test_deterministic_degree_values_counts_sum() {
        let m = 256u32;
        for dmax in 2..15 {
            let pdf = ideal_soliton_pdf(m as usize, dmax);
            for n in [0usize, 1, 7, 100, 1000] {
                let mut g = DeterministicDegreeValues::new_from_pdf(pdf.clone(), m);
                let seq = g.generate(n as u32);
                assert_eq!(seq.len(), n);
                let s: usize = g.counts().iter().sum();
                assert_eq!(s, n);
            }
        }
    }

    #[test]
    fn test_deterministic_additional_matches_full() {
        let pdf = vec![0u32, 90, 10];
        let m = 100u32;
        let mut g = DeterministicDegreeValues::new_from_pdf(pdf.clone(), m);
        let first = g.generate(10);
        let add = g.additional_degrees(10, 20);
        let mut comb = first.clone();
        comb.extend_from_slice(&add);
        comb.sort_unstable();

        let mut g2 = DeterministicDegreeValues::new_from_pdf(pdf, m);
        let full = g2.generate(20);
        assert_eq!(comb, full);
    }

    #[test]
    fn test_deterministic_from_cdf_ideal_soliton() {
        let m = 64usize;
        let dmax = 8usize;
        let cdf = ideal_soliton_cdf(m, dmax);
        let mut g = DeterministicDegreeValues::new_from_cdf(cdf);
        let seq = g.generate(50);
        assert_eq!(seq.len(), 50);
        assert_eq!(g.counts().iter().sum::<usize>(), 50);
    }
}