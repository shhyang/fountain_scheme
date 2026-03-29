// Copyright (c) 2025 Shenghao Yang. All rights reserved.
// Licensed under the MIT License. See LICENSE-MIT for details.

//! Random degree-set sampling in the spirit of RaptorQ (RFC 6330): separate **LT (active)** and
//! **PI (inactive)** index sets, using a CDF for LT degree and equal-distance sampling with a
//! prime modulus (see [`sample_degree_set_equal_distance`]).
//!
//! This module does **not** depend on the `rand` crate: randomness is supplied by [`PseudoRandom`].

use super::degree_dist::validate_cdf;
use super::degree_sampling::sample_degree_from_cdf;
use super::set_sampling::sample_degree_set_equal_distance;
use crate::math_utils::smallest_prime_ge;
use crate::types::PseudoRandom;
use fountain_engine::types::CodeParams;


/// RaptorQ-style generator: **W** LT symbols (active), **P** PI symbols (inactive), **L = W + P**.
///
/// - LT indices are in `[0, W)`; modulus uses `W' = smallest prime ≥ W` for equal-distance sampling.
/// - PI indices are in `[0, P)`; RFC 6330 uses `P1 = smallest prime ≥ L` for the PI walk — we use the
///   same `P1` here for [`sample_degree_set_equal_distance`].
/// - Inactive **degree** is `d1 ∈ {2, 3}` (matches RFC 6330 tuple constraints for the PI part).
pub struct RaptorDegreeSetGenerator<R: PseudoRandom> {
    cdf: Vec<u32>,
    /// Number of LT (active) intermediate symbols — `W` in RFC 6330.
    params: CodeParams,
    /// Smallest prime ≥ `W` (equal-distance modulus for active set).
    w_prime: usize,
    /// Smallest prime ≥ `P` (equal-distance modulus for inactive PI set).
    p_prime: usize,
    rng: R,
    is_systematic: bool,
}

impl<R: PseudoRandom> RaptorDegreeSetGenerator<R> {
    /// `cdf`: LT degree distribution (RFC Table 1 style). `w`, `p`: LT and PI symbol counts.
    pub fn new(cdf: Vec<u32>, params: CodeParams, rng: R) -> Self {
        assert!(!cdf.is_empty(), "cdf must be non-empty");
        let m = *cdf.last().expect("cdf non-empty");
        assert!(m != 0, "cdf last value (scale m) must be non-zero");
        if !validate_cdf(&cdf, m as usize) {
            panic!("Invalid CDF: {:?}", cdf);
        }
        let w_prime = smallest_prime_ge(params.num_active());
        let p_prime = smallest_prime_ge(params.num_pre_inactive());
        Self {
            cdf,
            params,
            w_prime,
            p_prime,
            rng,
            is_systematic: false,
        }
    }

    pub fn new_systematic(cdf: Vec<u32>, params: CodeParams, rng: R) -> Self {
        let mut generator = Self::new(cdf, params, rng);
        generator.is_systematic = true;
        generator
    }

    /// Sample `(active_lt_indices, inactive_pi_indices)` for one encoding symbol.
    ///
    /// Active indices are in `[0, W)`; inactive PI indices are in `[0, P)` (local PI indices; add `W`
    /// for global intermediate indices if your encoder uses a single buffer).
    pub fn degree_set(&mut self, coded_id: usize) -> Vec<usize> {
        let seed = u64::try_from(coded_id).expect("coded_id does not fit into u64");
        self.rng.seed(seed);
        let w = self.params.num_active();
        let p = self.params.num_pre_inactive();

        let m = *self.cdf.last().expect("cdf non-empty");
        let m_usize = usize::try_from(m).expect("CDF scale m must fit in usize");
        let r = self.rng.next(m_usize);
        let d = sample_degree_from_cdf(&self.cdf, r);

        let mut deg_set = if d == 0 {
            vec![]
        } else if w == 1 {
            vec![0]
        } else if self.is_systematic && coded_id < self.params.a {
            let d_lt = d.min(coded_id+1);
            if d_lt == 1 {
                vec![coded_id]
            } else if coded_id == 1 {
                vec![0, 1]
            } else {
                let w_sys_prime = smallest_prime_ge(coded_id);
                let interval_a = 1 + self.rng.next(w_sys_prime - 1);
                let start_a = self.rng.next(w_sys_prime);
                let mut set = sample_degree_set_equal_distance(w_sys_prime, coded_id, d_lt-1, start_a, interval_a);
                set.push(coded_id);
                set.sort();
                set
            }
        } else {
            let d_lt = d.min(w);
            let interval_a = 1 + self.rng.next(self.w_prime - 1);
            let start_a = self.rng.next(self.w_prime);
            sample_degree_set_equal_distance(self.w_prime, w, d_lt, start_a, interval_a)
        };

        let inactive_set = if p == 0 {
            Vec::new()
        } else {
            let d1 = (if d < 4 { 2 + self.rng.next(2) } else {2}).min(p);
            let interval_b = 1 + self.rng.next(self.p_prime - 1);
            let start_b = self.rng.next(self.p_prime);
            sample_degree_set_equal_distance(self.p_prime, p, d1, start_b, interval_b)
        };
            
        deg_set.extend(inactive_set.iter().map(|&i| i + w));
        deg_set
    }

    pub fn rng_mut(&mut self) -> &mut R {
        &mut self.rng
    }
}

pub struct R10DegreeSetGenerator<R: PseudoRandom> {
    cdf: Vec<u32>,
    k: usize,
    k_prime: usize,
    rng: R,
}

impl<R: PseudoRandom> R10DegreeSetGenerator<R> {
    pub fn new(cdf: Vec<u32>, k: usize, rng: R) -> Self {
        assert!(!cdf.is_empty(), "cdf must be non-empty");
        let m = *cdf.last().expect("cdf non-empty");
        assert!(m != 0, "cdf last value (scale m) must be non-zero");
        assert!(k > 0, "k must be positive");
        if !validate_cdf(&cdf, m as usize) {
            panic!("Invalid CDF: {:?}", cdf);
        }
        let k_prime = smallest_prime_ge(k);
        Self { cdf, k, k_prime, rng }
    }

    /// Sample active degree set for one coded symbol.
    ///
    /// - degree: `sample_degree_from_cdf(cdf, r)`
    /// - set: equal-distance sampling with random `(a, b)`
    pub fn degree_set(&mut self, coded_id: usize) -> Vec<usize> {
        let seed = u64::try_from(coded_id).expect("coded_id does not fit into u64");
        self.rng.seed(seed);

        let m = *self.cdf.last().expect("cdf non-empty");
        let m_usize = usize::try_from(m).expect("CDF scale m must fit in usize");
        let r = self.rng.next(m_usize);
        let d = sample_degree_from_cdf(&self.cdf, r);

        if self.k == 1 {
            return if d == 0 { vec![] } else { vec![0] };
        }
        if d == 0 {
            return vec![];
        }

        let d = d.min(self.k);
        let interval = 1 + self.rng.next(self.k_prime - 1);
        let start = self.rng.next(self.k_prime);
        sample_degree_set_equal_distance(self.k_prime, self.k, d, start, interval)
    }

    pub fn rng_mut(&mut self) -> &mut R {
        &mut self.rng
    }
}


#[cfg(test)]
#[allow(clippy::many_single_char_names)]
mod tests {
    use super::*;
    use crate::degree_sets::degree_sampling::sample_degree_from_cdf;
    use crate::validation::pseudo_rand::XorShift64;

    fn tiny_cdf() -> Vec<u32> {
        vec![0u32, 2, 4]
    }

    #[test]
    fn active_set_respects_k_and_bounds() {
        let cdf = tiny_cdf();
        let mut sampler = R10DegreeSetGenerator::new(cdf, 10, XorShift64::new(0x00C0_FFEE));
        for id in 0..50 {
            let set = sampler.degree_set(id);
            assert!(set.len() <= 10);
            assert!(set.iter().all(|&i| i < 10));
        }
    }

    #[test]
    fn rq_split_active_inactive_rfc6330_style_bounds() {
        let cdf = tiny_cdf();
        let k = 10;
        let a = 8;
        let l = 2;
        let h = 1;
        let params = CodeParams::new(k, a, l, h);
        let w = params.num_active();
        let p_pi = params.num_pre_inactive();
        let cap = w + p_pi;
        let mut g = RaptorDegreeSetGenerator::new(cdf, params, XorShift64::new(0xBEEF));
        for id in 0..30 {
            let set = g.degree_set(id);
            assert!(set.len() <= w + p_pi);
            assert!(
                set.iter().all(|&i| i < cap),
                "id {}: indices must lie in active [0,{w}) or PI [w,{cap})",
                id
            );
        }
    }

    // check that RaptorDegreeSetGenerator generates the same output as R10DegreeSetGenerator
    #[test]
    fn systematic_triangular_degree_matches_d_eff() {
        let cdf = tiny_cdf();
        // `p = 0` so `degree_set` has no `+ w` inactive tail; we only check the active triangular part.
        let k = 10;
        let a = 10;
        let l = 0;
        let h = 0;
        let w = a + l;
        let mut g =
            RaptorDegreeSetGenerator::new_systematic(cdf.clone(), CodeParams::new(k, a, l, h), XorShift64::new(0xACE));
        for coded_id in 0..a {
            let set = g.degree_set(coded_id);
            assert!(set.iter().all(|&i| i < w), "id {} set {:?} out of range", coded_id, set);
            if coded_id == 0 {
                assert_eq!(set, vec![0]);
                continue;
            }
            if coded_id == 1 {
                let mut r = XorShift64::new(0xACE);
                r.seed(1);
                let m = *cdf.last().unwrap() as usize;
                let sample_d = sample_degree_from_cdf(&cdf, r.next(m));
                let d_lt = sample_d.min(coded_id + 1);
                if d_lt == 1 {
                    assert_eq!(set, vec![1]);
                } else {
                    assert_eq!(set, vec![0, 1]);
                }
                continue;
            }
            let mut r = XorShift64::new(0xACE);
            r.seed(coded_id as u64);
            let m = *cdf.last().unwrap() as usize;
            let sample_d = sample_degree_from_cdf(&cdf, r.next(m));
            let d_eff = sample_d.min(coded_id + 1);
            if d_eff == 0 {
                assert!(set.is_empty(), "d>=1 and coded_id>=2 implies d_eff>=1");
            } else if d_eff == 1 {
                assert_eq!(set, vec![coded_id]);
            } else {
                assert_eq!(set.len(), d_eff, "id {} set {:?}", coded_id, set);
                assert_eq!(*set.last().unwrap(), coded_id);
                assert!(set[..set.len() - 1].iter().all(|&i| i < coded_id));
            }
        }
    }

    #[test]
    fn cross_verify_r10_rq() {
        let cdf = tiny_cdf();
        let k = 10;
        let a = 8;
        let l = 2;
        let h = 1;
        let w = a + l;
        let mut g = RaptorDegreeSetGenerator::new(cdf.clone(), CodeParams::new(k, a, l, h), XorShift64::new(0xBEEF));
        let mut g1 = R10DegreeSetGenerator::new(cdf, w, XorShift64::new(0xBEEF));
        for id in 0..30 {
            let merged = g.degree_set(id);
            let active_lt: Vec<usize> = merged.iter().copied().filter(|&i| i < w).collect();
            let active1 = g1.degree_set(id);
            assert_eq!(active_lt, active1, "LT active part should match R10 on same W");
        }
    }
}
