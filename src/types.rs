// Copyright (c) 2025 Shenghao Yang. All rights reserved.
// Licensed under the MIT License. See LICENSE-MIT for details.


/// New types for probability values and degree values 
pub type ProbValue = u32;

/// Generic pseudo-random source without depending on `rand`.
///
/// - `seed` uses `u64` so sequences are stable across 32/64-bit targets.
/// - `next(hi)` returns a value uniformly distributed in `[0, hi)` (requires `hi > 0`).
pub trait PseudoRandom: Clone {
    fn seed(&mut self, seed: u64);
    /// Uniform integer in `[0, hi)`; `hi` must be positive.
    fn next(&mut self, hi: usize) -> usize;
}