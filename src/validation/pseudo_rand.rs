// Copyright (c) 2026 Shenghao Yang.
// All rights reserved.

use crate::types::PseudoRandom;

/// Small dependency-free PRNG for tests and examples (xorshift64).
#[derive(Clone, Debug)]
pub struct XorShift64 {
    state: u64,
}

impl XorShift64 {
    #[must_use]
    pub fn new(seed: u64) -> Self {
        Self { state: seed.max(1) }
    }

    fn next_u64(&mut self) -> u64 {
        let mut x = self.state;
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        self.state = x;
        x
    }
}
impl PseudoRandom for XorShift64 {
    fn seed(&mut self, seed: u64) {
        self.state = seed.max(1);
    }

    fn next(&mut self, hi: usize) -> usize {
        assert!(hi > 0, "next: hi must be positive");
        let x = self.next_u64();
        usize::try_from((u128::from(x) * hi as u128) >> 64).expect("hi too large for usize")
    }
}
