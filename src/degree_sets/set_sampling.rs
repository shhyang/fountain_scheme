// Copyright (c) 2025 Shenghao Yang. All rights reserved.
// Licensed under the MIT License. See LICENSE-MIT for details.

/// Sample the degree set with a sequence of sampled values.
/// Given a sequence of \(d\) values, sample the degree set with the values.
/// The \(d\) sampled values have the property that sample[i] is sampled from [0, n-1-i).
/// Returns a vector of \(d\) values.
pub fn sample_degree_set_with_values(k: usize, sample: &[usize]) -> Vec<usize> {
    let d = sample.len();
    if d > k {
        panic!("sample length must be less than or equal to n");
    }
    let mut values: Vec<usize> = (0..k).collect();
    for (step, &r) in sample.iter().enumerate() {
        let remaining = k - step;
        if r >= remaining {
            panic!("sample[{step}] must be in 0..{remaining}, got {r}");
        }
        values.swap(r, remaining - 1);
    }
    let mut set = values[k - d..].to_vec();
    set.sort();
    set
}

/// Equal-distance sampling of degree sets.
/// Given a sequence of \(n\) values, sample \(d\) values with equal distance.
/// The values are sampled from \(a, a+b, a+2b, \ldots, a+(d-1)b\), 
/// where \(a\) is uniformly distributed in \(0,1,\dots, n-1\) and \(b\) is uniformly distributed in \(1,2,\dots, n-1\).
/// Returns a vector of \(d\) values.
/// To avoid duplications, \(d\) should not divide \(n\).
pub fn sample_degree_set_equal_distance(k_prime: usize, k: usize, d: usize, start: usize, interval: usize) -> Vec<usize> {
    // valid the input
    if d > k {
        panic!("d must be less than or equal to n");
    }
    if d == 0 {
        return vec![];
    }
    if interval == 0 {
        panic!("interval must be greater than 0");
    }
    let mut start = start;
    while start >= k {
        start = (start + interval) % k_prime;
    }
    let mut set = Vec::with_capacity(d);
    set.push(start);
    for _ in 1..d {
        start = (start + interval) % k_prime;
        while start >= k {
            start = (start + interval) % k_prime;
        }
        set.push(start);
    }
    set.sort();
    set
}



/// Cyclic-stride approach to sample the degree set without random number generator.
pub struct CyclicStrideDegreeSets {
    k: usize,
    a: usize,
    b: usize,
    p: usize,
}

impl CyclicStrideDegreeSets {
    /// Create a cyclic-stride deterministic sampler.
    ///
    /// The generated sequence is `(a + t*b) mod k` with a moving offset.
    /// To avoid duplication within each sampled set (when `d <= k`), require `gcd(b, k) = 1`.
    pub fn new(k: usize) -> Self {
        if k == 0 {
            panic!("k must be greater than 0");
        }
        Self { k, a: 0, b: 1, p: 0 }
    }

    /// Sample one degree set with size `d` in `O(d)`.
    pub fn sample_set(&mut self, d: usize) -> Vec<usize> {
        if d > self.k {
            panic!("d must be less than or equal to k");
        }
        let mut set = Vec::with_capacity(d);
        for t in 0..d {
            let q = (self.p + t) % self.k;
            let idx =
                ((self.a as u128 + (q as u128 * self.b as u128)) % self.k as u128) as usize;
            set.push(idx);
        }
        self.p += d;
        while self.p >= self.k {
            self.a = (self.a + 1) % self.k;
            self.b = (self.b + 1) % self.k;
            while gcd_usize(self.b, self.k) != 1 {
                self.b = (self.b + 1) % self.k;
            }
            self.p -= self.k;
        }
        set.sort();
        set
    }

    /// Sample multiple degree sets from a degree sequence.
    pub fn sample_sets(&mut self, degrees: &[usize]) -> Vec<Vec<usize>> {
        degrees.iter().map(|&d| self.sample_set(d)).collect()
    }
}

fn gcd_usize(mut a: usize, mut b: usize) -> usize {
    while b != 0 {
        let r = a % b;
        a = b;
        b = r;
    }
    a
}

