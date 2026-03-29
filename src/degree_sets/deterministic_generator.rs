// Copyright (c) 2025 Shenghao Yang. All rights reserved.
// Licensed under the MIT License. See LICENSE-MIT for details.

use super::degree_sampling::DeterministicDegreeValues;
use super::set_sampling::CyclicStrideDegreeSets;

/// Deterministic degree-set generator: [`DeterministicDegreeValues`] for per-symbol degrees and
/// [`CyclicStrideDegreeSets`] for index sets.
///
/// **Caching:** Each block of `k` consecutive coded indices loads `k` degree sets at once via
/// [`CyclicStrideDegreeSets::sample_sets`], so the internal stride pointer matches emitting those
/// symbols in order.
///
/// **Contract:** For the canonical global stride schedule, call [`Self::degree_set`] with
/// `coded_id = 0, 1, 2, …` in order. The first call loads the block `[0, k)`; if the first call
/// uses another `coded_id`, earlier indices are skipped and the stride state does not include them.
pub struct CyclicStrideDegreeSetGenerator {
    k: usize, 
    degree_values: DeterministicDegreeValues,
    sampler: CyclicStrideDegreeSets,
    current_range: (usize, usize),
    current_degree_set: Vec<Vec<usize>>,
}

impl CyclicStrideDegreeSetGenerator {
    /// `k` is both the index universe size (valid indices are `0..k`) and the block length for degree batches.
    /// Starts with an empty cached block; the first [`Self::degree_set`] call loads the block
    /// that contains `coded_id` (typically start with `coded_id == 0`).
    pub fn new(k: usize, cdf: Vec<u32>) -> Self {
        let degree_values = DeterministicDegreeValues::new_from_cdf(cdf);
        let sampler = CyclicStrideDegreeSets::new(k);
        Self { k, degree_values, sampler, current_range: (0, 0), current_degree_set: Vec::new() }
    }

    /// Returns the active index set for `coded_id`, cloning from the current block cache when hit.
    pub fn degree_set(&mut self, coded_id: usize) -> Vec<usize> {
        let (start, end) = self.current_range;
        if coded_id >= start && coded_id < end {
            self.current_degree_set[coded_id - start].clone()
        } else {
            let degree = self.degree_values.additional_degrees(coded_id as u32, coded_id as u32 + self.k as u32);
            self.current_range = (coded_id, coded_id + self.k);
            self.current_degree_set = self.sampler.sample_sets(&degree);
            self.current_degree_set[0].clone()
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Minimal valid CDF: degrees 1..=2 with mass 4 (see degree_sampling tests).
    fn tiny_cdf() -> Vec<u32> {
        vec![0, 2, 4]
    }

    #[test]
    fn cyclic_stride_generator_first_block_cached() {
        let k = 4;
        let mut g = CyclicStrideDegreeSetGenerator::new(k, tiny_cdf());
        for id in 0..k {
            let s = g.degree_set(id);
            assert!(
                s.len() <= k,
                "degree set size should not exceed k"
            );
        }
    }

    #[test]
    fn cyclic_stride_generator_second_block() {
        let k = 4;
        let mut g = CyclicStrideDegreeSetGenerator::new(k, tiny_cdf());
        for id in 0..k {
            let _ = g.degree_set(id);
        }
        let _ = g.degree_set(k);
    }
}