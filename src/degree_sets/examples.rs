// Copyright (c) 2025 Shenghao Yang. All rights reserved.
// Licensed under the MIT License. See LICENSE-MIT for details.

use super::traits::DegreeSetGenerator;
use fountain_engine::CodeParams;

/// Example custom degree set that generates a simple deterministic pattern
pub struct DeterministicDegreeSet {
    params: CodeParams,
}

impl DeterministicDegreeSet {
    /// Creates a new deterministic degree set from the given code parameters.
    pub fn new(params: &CodeParams) -> Self {
        Self {
            params: params.clone(),
        }
    }
}

impl DegreeSetGenerator for DeterministicDegreeSet {
    fn degree_set(&mut self, coded_id: usize) -> (Vec<usize>, Vec<usize>) {
        // Simple deterministic pattern: each coded vector connects to exactly one message vector
        let active_indices = if coded_id < self.params.k {
            vec![coded_id]
        } else {
            vec![coded_id % self.params.k]
        };
        (active_indices, vec![])
    }
    fn name(&self) -> &str {
        "Deterministic"
    }
}

/// Example custom degree set that generates a fixed degree pattern
pub struct FixedDegreeSet {
    params: CodeParams,
    degree: usize,
}

impl FixedDegreeSet {
    /// Creates a new fixed-degree set with the given degree and code parameters.
    pub fn new(params: &CodeParams, degree: usize) -> Self {
        Self {
            params: params.clone(),
            degree,
        }
    }
}

impl DegreeSetGenerator for FixedDegreeSet {
    fn degree_set(&mut self, coded_id: usize) -> (Vec<usize>, Vec<usize>) {
        use rand::seq::SliceRandom;
        use rand::{SeedableRng, rngs::StdRng};

        let mut rng = StdRng::seed_from_u64(coded_id as u64);
        let mut indices: Vec<usize> = (0..self.params.k).collect();
        indices.shuffle(&mut rng);
        indices.truncate(self.degree);
        indices.sort();

        (indices, vec![])
    }
    fn name(&self) -> &str {
        "FixedDegree"
    }
}

/// Deterministic triangular degree set: coded vector `i` connects to source symbols `0..=i`.
pub struct CustomTriangularDegreeSet {
    params: CodeParams,
}

impl CustomTriangularDegreeSet {
    /// Creates a new custom triangular degree set from the given code parameters.
    pub fn new(params: &CodeParams) -> Self {
        Self {
            params: params.clone(),
        }
    }
}

impl DegreeSetGenerator for CustomTriangularDegreeSet {
    fn degree_set(&mut self, coded_id: usize) -> (Vec<usize>, Vec<usize>) {
        if coded_id == 0 {
            return (vec![0], vec![]);
        }

        // Triangular pattern: degree increases with coded_id
        let degree = (coded_id % self.params.k).min(self.params.k - 1) + 1;
        let mut indices = Vec::new();

        for i in 0..degree {
            indices.push(i);
        }

        (indices, vec![])
    }
    fn name(&self) -> &str {
        "CustomTriangular"
    }
}
