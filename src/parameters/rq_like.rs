// Copyright (c) 2025 Shenghao Yang. All rights reserved.
// Licensed under the MIT License. See LICENSE-MIT for details.

use crate::math_utils::{is_prime, smallest_prime_ge};
use crate::parameters::r10::smallest_x_for_ldpc;
use fountain_engine::types::CodeParams;

/// Generate "RaptorQ-like" parameters for arbitrary `k`.
///
/// This is an approximation rule for cases where an exact RFC table lookup is not used:
/// - `l` is chosen as a prime number near `0.01k + sqrt(2k)`.
/// - `h` is chosen in `[10, 16]`.
/// - `a` is chosen so that `a + l` is prime and `k - a + h` is proportional to `sqrt(k)`.
///
/// Note: RFC6330 uses a fixed table and next-larger lookup for exact interoperability.
pub fn generate_rq_like_parameters(k: usize) -> CodeParams {
    assert!(k >= 1, "k must be >= 1");

    let x = smallest_x_for_ldpc(k);
    let l_target = (0.01 * k as f64 + x as f64).ceil() as usize;
    let l = smallest_prime_ge(l_target.max(2));

    // Keep h in [10,16], increasing slowly with k.
    let h = choose_h_10_to_16(k);

    // Piecewise rule from docs:
    // - k < 60: a is almost k
    // - k >= 60: a ~= k - 1.57*sqrt(k)
    let target_a = if k < 60 {
        k
    } else {
        let reduction = (1.57 * (k as f64).sqrt()).round() as usize;
        k.saturating_sub(reduction).max(1)
    };
    let a = nearest_a_with_prime_sum(target_a, k, l);

    CodeParams::new(k, a, l, h)
}

fn choose_h_10_to_16(k: usize) -> usize {
    // Map k in [1, 56403] approximately to [10,16]. For larger k, clamp to 16.
    let min_k: f64 = 60.0_f64.ln();
    let max_k: f64 = 56_403.0_f64.ln();
    let h = 10.0 + 6.0 * ((k as f64).ln() - min_k) / (max_k - min_k);
    h.round().clamp(10.0, 16.0) as usize
}

fn nearest_a_with_prime_sum(target_a: usize, k: usize, l: usize) -> usize {
    let t = target_a.clamp(1, k);
    if is_prime(t + l) {
        return t;
    }
    // Search outward from target and return the closest valid `a`.
    for d in 1..=k {
        if t + d <= k {
            let a_up = t + d;
            if is_prime(a_up + l) {
                return a_up;
            }
        }
        if t > d {
            let a_down = t - d;
            if is_prime(a_down + l) {
                return a_down;
            }
        }
    }
    // Fallback (should not happen in practice for these ranges).
    1
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_rq_like_basic_properties() {
        // k == 1 is degenerate: only a = 1 is possible, but with prime l ≥ 3 the sum 1 + l is
        // even and > 2, so this heuristic cannot satisfy is_prime(a + l).
        for k in [2usize, 10, 100, 1000, 10_000, 56_403] {
            let p = generate_rq_like_parameters(k);
            assert_eq!(p.k, k);
            assert!((10..=16).contains(&p.h));
            assert!(is_prime(p.l));
            assert!(is_prime(p.a + p.l));
            assert!(p.a >= 1 && p.a <= k);
        }
    }

    #[test]
    fn test_rq_like_l_close_to_target() {
        for k in [100usize, 1000, 10_000, 56_403] {
            let p = generate_rq_like_parameters(k);
            let target = 0.01 * k as f64 + (2.0 * k as f64).sqrt();
            let diff = (p.l as f64 - target).abs();
            assert!(diff < 20.0);
        }
    }

    #[test]
    fn test_rq_like_a_piecewise_rule() {
        for k in [10usize, 30, 59] {
            let p = generate_rq_like_parameters(k);
            assert!(k.saturating_sub(p.a) <= 5, "k<60 should keep a near k");
        }
        for k in [60usize, 100, 1000, 10_000] {
            let p = generate_rq_like_parameters(k);
            let target = k as f64 - 1.57 * (k as f64).sqrt();
            assert!((p.a as f64 - target).abs() <= 10.0);
        }
    }

    // print the parameters for a range of k
    #[test]
    fn test_print_parameters_for_range_of_k() {
        for k in [10, 30, 49, 60, 101, 180, 168, 573, 811, 1136, 1579, 2125, 2831, 3751, 4901, 6296, 8111, 10351, 13143, 16674, 21199, 26838, 33961, 42916, 56403] {
            let p = generate_rq_like_parameters(k);
            println!("k: {}, a: {}, l: {}, h: {}", k, p.a, p.l, p.h);
        }
    }
}




