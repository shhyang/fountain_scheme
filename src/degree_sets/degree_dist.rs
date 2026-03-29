// Copyright (c) 2025 Shenghao Yang. All rights reserved.
// Licensed under the MIT License. See LICENSE-MIT for details.

//! Soliton-family degree distributions used by fountain codes.
//!
//! This module supports two representations of a degree distribution over degrees \(d=1..D_{max}\):
//!
//! - **Floating point PDF**: `Vec<f64>` where `pdf[d]` is \(p_d\) and `pdf[0]=0`.
//! - **Fixed-point CDF/PDF**: probabilities represented as integers with an `accuracy` in bits.
//!   Let \(A = 2^{accuracy}-1\). We represent:
//!   - `pdf[d]` as an integer in \([0, A]\) with \(\sum_d pdf[d] = A\).
//!   - `cdf[d] = \sum_{i=1}^d pdf[i]\), with `cdf[0]=0` and `cdf[dmax]=A`.
//!
//! This matches the common convention \(p_d = c_d - c_{d-1}\) with \(c_0=0\), and provides a
//! deterministic way to sample degrees using integer arithmetic.

use crate::types::ProbValue;

pub fn ideal_soliton_pdf(m: usize, dmax: usize) -> Vec<ProbValue> {
    let cdf = ideal_soliton_cdf(m, dmax);
    cdf_to_pdf(&cdf)
}

pub fn robust_soliton_pdf(m: usize, dmax: usize, c: f64, delta: f64) -> Vec<ProbValue> {
    let cdf = robust_soliton_cdf(m, dmax, c, delta);
    cdf_to_pdf(&cdf)
}

pub fn ideal_soliton_cdf(m: usize, dmax: usize) -> Vec<ProbValue> {
    let mut cdf = vec![0usize; dmax+1];
    cdf[1] = m / dmax;
    for i in 2..dmax {
        cdf[i] = m - m / i + cdf[1];
    }
    cdf[dmax] = m;
    cdf.iter().map(|x| *x as u32).collect()
}

pub fn robust_soliton_cdf(m: usize, dmax: usize, c: f64, delta: f64) -> Vec<ProbValue> {
    let ideal_pdf = ideal_soliton_pdf(m, dmax)
        .iter()
        .map(|x| *x as f64 / m as f64)
        .collect::<Vec<f64>>();
    
    let kf = dmax as f64;
    let r = (kf.sqrt() * c * (kf / delta).ln()).max(1.0);
    let upper = (kf / r).ceil() as usize;

    let mut tau = vec![0.0; dmax + 1];
    for i in 1..upper.min(dmax) {
        tau[i] = r / (i as f64 * kf);
    }
    if upper <= dmax {
        tau[upper] = r * (r / delta).ln() / kf;
    }

    let mut cdf = vec![0.0; dmax + 1];
    for i in 1..=dmax {
        cdf[i] = cdf[i-1] + ideal_pdf[i] + tau[i];
    }
    let a = cdf[dmax] as f64;
    let mut cdf32 = vec![0u32; dmax + 1];
    for i in 1..dmax {
        cdf32[i] = (cdf[i] / a * m as f64) as u32;
    }
    cdf32[dmax] = m as u32;
    if cdf32[dmax-1] > cdf32[dmax] {
        panic!("Invalid degree distribution: cdf[dmax-1] > cdf[dmax]");
    }
    cdf32
}

/// Generate the CDF of the asymptotically optimal degree distribution.
///
/// `m` is the total probability mass (same scale as [`validate_cdf`]: last entry must be `m`).
/// `d1` is the cumulative mass through degree 1; the remaining mass `m - d1` follows the
/// ideal-soliton-style tail `m2 - m2/i` for degrees `2..dmax-1`, then `cdf[dmax] = m`.
pub fn asymp_optimal_cdf(m: usize, d1: usize, dmax: usize) -> Vec<ProbValue> {
    if d1 > m {
        panic!("Invalid degree distribution: d1 > m");
    }
    if dmax < 1 {
        panic!("dmax must be at least 1");
    }
    let mut cdf = vec![0usize; dmax + 1];
    cdf[1] = d1;
    let m2 = m - d1;
    for i in 2..dmax {
        cdf[i] = d1 + m2 - m2 / i;
    }
    cdf[dmax] = m;
    cdf.iter().map(|x| *x as u32).collect()
}

/// Implement the Raptor10 and RaptorQ degree distribution.
pub const RAPTOR10_CDF: [ProbValue; 8] = [
    0, 10241, 491582, 712794, 831695, 948446, 1032189, 1048576
];

pub const RAPTORQ_CDF: [ProbValue; 31] = [
    0, 5243, 529531, 704294, 791675, 844104, 879057, 904023, 
    922747, 937111, 948962, 958494, 966438, 973160, 978921,
    983914, 988283, 992138, 995565, 998631, 1001391, 1003887, 
    1006157, 1008229, 1010129, 1011876, 1013490, 1014983, 1016370, 
    1017662, 1048576
];

pub fn cdf_to_pdf(cdf: &[ProbValue]) -> Vec<ProbValue> {
    let mut pdf = vec![0u32; cdf.len()];
    for i in 1..cdf.len() {
        pdf[i] = cdf[i] - cdf[i-1];
    }
    pdf
}

pub fn pdf_to_cdf(pdf: &[ProbValue]) -> Vec<ProbValue> {
    let mut cdf = vec![0u32; pdf.len()];
    for i in 1..pdf.len() {
        cdf[i] = cdf[i-1] + pdf[i];
    }
    cdf
}

/// Validate the PDF of a degree distribution.
/// The PDF is a vector of integers of `m` bits, 
/// where the sum of the PDF is equal to the largest integer of `m` bits.
/// For a degree distribution, we always have `pdf[0] = 0`.
pub fn validate_pdf(pdf: &[ProbValue], m: usize) -> bool {
    if pdf[0] != 0 {
        return false;
    }
    let sum: u64 = pdf.iter().map(|&x| x as u64).sum();
    if sum == 0 {
        return false;
    }
    if sum != m as u64 {
        return false;
    }
    true
}

/// Validate the CDF of a degree distribution.
/// The CDF is a vector of integers of `m` bits, 
/// where the last element is the largest integer of `m` bits.
/// The first element is 0.
/// The CDF is non-decreasing.
pub fn validate_cdf(cdf: &[ProbValue], m: usize) -> bool {  
    if cdf.len() < 2 {
        return false;
    }
    if cdf[0] != 0 {
        return false;
    }
    if cdf[cdf.len() - 1] != m as u32 {
        return false;
    }
    if m == 0 {
        return false;
    }

    // Invariants: 0 <= c_d <= c_{d+1} <= max_value
    for win in cdf.windows(2) {
        let a = win[0];
        let b = win[1];
        if a > m as u32 || b > m as u32 {
            return false;
        }
        if a > b {
            return false;
        }
    }

    true
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::rngs::StdRng;
    use rand::Rng;
    use rand::SeedableRng;

    #[test]
    fn test_raptor10_cdf() {
        assert!(validate_cdf(RAPTOR10_CDF.as_slice(), 1048576));
    }

    #[test]
    fn test_raptorq_cdf() {
        assert!(validate_cdf(RAPTORQ_CDF.as_slice(), 1048576));
    }

    #[test]
    fn print_cdf() {
        let m = 2usize << 19;
        //println!("m: {}", m);
        let cdf = asymp_optimal_cdf(m, 10241, 7);
        println!("cdf of ideal soliton distribution with dmax=7: {:?}", cdf);
        let cdf = asymp_optimal_cdf(m, 5243, 31);
        println!("cdf of ideal soliton distribution with dmax=31: {:?}", cdf);
    }

    #[test]
    fn test_asymp_optimal_cdf() {
        for m in 32..=1024 {
            for d1 in 0..=10 {
                let cdf = asymp_optimal_cdf(m, d1, 100);
                assert!(validate_cdf(&cdf, m));
            }
        }
    }

    /// Generate a random integer PDF with `entries` bins whose sum is `2^m - 1`.
    /// For degree distributions, `pdf[0]` is fixed to 0.
    fn random_pdf(entries: usize, m: usize, rng: &mut StdRng) -> Vec<u32> {
        assert!(entries >= 2, "entries must be at least 2");
        assert!((1..=32).contains(&m), "m must be in 1..=32");

        let mut pdf = vec![0u32; entries];
        let mut remaining = m as u32;

        // Randomly split `remaining` across entries 1..entries-1.
        for item in pdf.iter_mut().take(entries - 1).skip(1) {
            let v = rng.gen_range(0..=remaining);
            *item = v;
            remaining -= v;
        }
        pdf[entries - 1] = remaining;
        pdf
    }

    #[test]
    fn test_cdf_to_pdf() {
        let mut rng = StdRng::seed_from_u64(20_260_315);
        for m in [8usize, 12, 16, 24, 32] {
            for entries in [2usize, 3, 5, 10, 32, 64] {
                for _ in 0..20 {
                    let pdf = random_pdf(entries, m, &mut rng);
                    assert!(validate_pdf(&pdf, m));

                    let cdf = pdf_to_cdf(&pdf);
                    let new_pdf = cdf_to_pdf(&cdf);
                    let new_cdf = pdf_to_cdf(&new_pdf);

                    assert_eq!(new_cdf, cdf);
                    assert_eq!(new_pdf, pdf);
                }
            }
        }
    }

    #[test]
    fn test_ideal_soliton_pdf() {
        for m in 8..=32 {
            for k in 1..100 {
                let pdf = ideal_soliton_pdf(m, k);
                assert!(validate_pdf(&pdf, m));
            }
        }
    }
    #[test]
    fn test_robust_soliton_pdf() {
        for m in 8..=32 {
            for k in 1..100 {
                let pdf = robust_soliton_pdf(m, k, 0.1, 0.05);
                assert!(validate_pdf(&pdf, m));
            }
        }
    }

    #[test]
    fn test_ideal_soliton_cdf() {
        for m in 8..=32 {
            for dmax in 1..100 {
                let cdf = ideal_soliton_cdf(m, dmax);
                assert!(validate_cdf(&cdf, m));
            }
        }
    }

    #[test]
    fn test_robust_soliton_cdf() {
        for m in 8..=32 {
            for dmax in 1..100 {
                let cdf = robust_soliton_cdf(m, dmax, 0.1, 0.05);
                assert!(validate_cdf(&cdf, m));
            }
        }
    }

    #[test]
    fn test_validate_cdf_invalid() {
        let dmax = 10;
        let m = 1usize << 8;

        let base = ideal_soliton_cdf(m, dmax);

        // c0 must be 0 under the current validate_cdf implementation
        let mut bad = base.clone();
        bad[0] = 1;
        assert!(!validate_cdf(&bad, m));

        // must be non-decreasing
        let mut bad = base.clone();
        let prev = bad[1];
        bad[2] = prev.saturating_sub(1); // force a potential decrease
        assert!(!validate_cdf(&bad, m));

        // last value must be max_value
        let mut bad = base.clone();
        bad[dmax] = (m - 1) as u32;
        assert!(!validate_cdf(&bad, m));

        // entries must be <= max_value
        let mut bad = base.clone();
        bad[1] = (m + 1) as u32;
        assert!(!validate_cdf(&bad, m));
    }
}