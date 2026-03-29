// Copyright (c) 2025 Shenghao Yang. All rights reserved.
// Licensed under the MIT License. See LICENSE-MIT for details.

/// Compute the binomial coefficient C(n, r) in `u128`.
///
/// Uses the multiplicative formula and symmetry (`r = min(r, n-r)`).
pub fn binomial(n: usize, r: usize) -> u128 {
    if r > n {
        return 0;
    }
    let r = r.min(n - r);
    if r == 0 {
        return 1;
    }
    let mut result = 1u128;
    for i in 1..=r {
        let num = (n - r + i) as u128;
        let den = i as u128;
        result = (result * num) / den;
    }
    result
}

/// Compute the middle binomial coefficient C(h, ceil(h/2)).
pub fn middle_binomial(h: usize) -> u128 {
    let r = h.div_ceil(2);
    binomial(h, r)
}

/// Return the smallest prime number greater than or equal to `n`.
pub fn smallest_prime_ge(mut n: usize) -> usize {
    if n <= 2 {
        return 2;
    }
    if n.is_multiple_of(2) {
        n += 1;
    }
    while !is_prime(n) {
        n += 2;
    }
    n
}

/// Primality test for positive integers.
pub fn is_prime(n: usize) -> bool {
    if n < 2 {
        return false;
    }
    if n == 2 {
        return true;
    }
    if n.is_multiple_of(2) {
        return false;
    }
    let mut d = 3usize;
    while d * d <= n {
        if n.is_multiple_of(d) {
            return false;
        }
        d += 2;
    }
    true
}

