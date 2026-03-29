// Copyright (c) 2025 Shenghao Yang. All rights reserved.
// Licensed under the MIT License. See LICENSE-MIT for details.

use fountain_engine::types::CodeParams;
use crate::math_utils::{middle_binomial, smallest_prime_ge};

/// Generate Raptor10-like parameters from `k` source symbols.
///
/// Follows the rules documented in `docs/doc-scheme.org`:
/// - no pre-inactivation: `a = k`, `i = 0`
/// - `x` is the smallest positive integer with `x(x-1) >= 2k`
/// - `l` is the smallest prime number `>= 0.01k + x`
/// - `h` is the smallest integer satisfying `k + l <= C(h, ceil(h/2))`
pub fn generate_r10_parameters(k: usize) -> CodeParams {
    let x = smallest_x_for_ldpc(k);
    let l_target = (0.01 * k as f64 + x as f64).ceil() as usize;
    let l = smallest_prime_ge(l_target.max(2));
    let h = smallest_h_for_hdpc(k + l);
    CodeParams::new(k, k, l, h)
}

pub(crate) fn smallest_x_for_ldpc(k: usize) -> usize {
    let kf = k as f64;
    let mut x = ((1.0 + (1.0 + 8.0 * kf).sqrt()) / 2.0).ceil() as usize;
    if x == 0 {
        x = 1;
    }
    while x > 1 && (x - 1).saturating_mul(x.saturating_sub(2)) >= 2usize.saturating_mul(k) {
        x -= 1;
    }
    while x.saturating_mul(x.saturating_sub(1)) < 2usize.saturating_mul(k) {
        x += 1;
    }
    x
}

fn smallest_h_for_hdpc(target_cols: usize) -> usize {
    // Find the smallest h such that C(h, ceil(h/2)) >= target_cols.
    let mut h = 1usize;
    while middle_binomial(h) < target_cols as u128 {
        h += 1;
    }
    h
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::math_utils::{binomial, is_prime};

    #[test]
    fn test_r10_params_no_pre_inactivation() {
        for k in [1usize, 10, 100, 1000] {
            let p = generate_r10_parameters(k);
            assert_eq!(p.k, k);
            assert_eq!(p.a, k);
            assert_eq!(p.b, 0);
        }
    }

    #[test]
    fn test_r10_ldpc_rule() {
        for k in [1usize, 10, 55, 256, 1000] {
            let p = generate_r10_parameters(k);
            let x = smallest_x_for_ldpc(k);
            assert!(x * (x - 1) >= 2 * k);
            if x > 1 {
                assert!((x - 1) * (x - 2) < 2 * k);
            }
            assert!(is_prime(p.l));
            assert!(p.l as f64 >= 0.01 * k as f64 + x as f64);
        }
    }

    #[test]
    fn test_r10_hdpc_rule() {
        for k in [1usize, 10, 55, 256, 1000] {
            let p = generate_r10_parameters(k);
            let h_prime = p.h.div_ceil(2);
            assert!(binomial(p.h, h_prime) >= (k + p.l) as u128);
            if p.h > 1 {
                let prev_h = p.h - 1;
                let prev_h_prime = prev_h.div_ceil(2);
                assert!(binomial(prev_h, prev_h_prime) < (k + p.l) as u128);
            }
        }
    }
}
