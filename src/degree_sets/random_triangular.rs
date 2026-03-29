// Copyright (c) 2025 Shenghao Yang. All rights reserved.
// Licensed under the MIT License. See LICENSE-MIT for details.

use super::SolitonDistribution;
use super::sampling::*;
use super::traits::DegreeSetGenerator;
use fountain_engine::CodeParams;
use rand::{SeedableRng, rngs::StdRng};

/// Degree set using random triangular embedding.
pub struct RandomTriangularDegreeSet {
    active_sampling: DegreeSampling,
    inactive_sampling: DegreeSampling,
    //params: CodeParams, // range of indices for sampling
    k: usize,
    num_active: usize,
    num_inactive: usize,
}

impl RandomTriangularDegreeSet {
    /// Creates a new random-triangular degree set from the distribution and code parameters.
    pub fn new(degree_distribution: &SolitonDistribution, params: &CodeParams) -> Self {
        let mut active_pdf = degree_distribution.pdf();
        shorten_pdf(&mut active_pdf, params.num_active());
        let active_cdf = pdf_to_cdf(&active_pdf);

        // the degree of inactive packets is 2 or 3
        let mut inactive_pdf = vec![0.0, 0.0, 0.5, 0.5];
        shorten_pdf(&mut inactive_pdf, params.num_inactive());
        let inactive_cdf = pdf_to_cdf(&inactive_pdf);

        let active_sampling = DegreeSampling::new(active_cdf);
        let inactive_sampling = DegreeSampling::new(inactive_cdf);
        Self {
            active_sampling,
            inactive_sampling,
            k: params.k,
            num_active: params.num_active(),
            num_inactive: params.num_inactive(),
        }
    }

    /// Generates a random lower-triangular index set for the given symbol index.
    ///
    /// Like [`super::triangular::TriangularDegreeSet`], the returned set always includes `idx` as its maximum,
    /// but the remaining indices are chosen uniformly at random.
    pub fn triangular_random(&mut self, seeded_rng: &mut StdRng, idx: usize) -> Vec<usize> {
        if idx == 0 {
            return vec![0];
        }
        if idx == 1 {
            return vec![0, 1];
        }
        // for coded_id = 2,...,k_a-1
        let degree = self.active_sampling.random(seeded_rng).min(idx);

        if degree == 0 {
            return vec![];
        }
        if degree == 1 {
            return vec![idx];
        }
        let mut indices = sampling_indices_with_degree(degree - 1, idx - 1, seeded_rng);
        indices.push(idx);
        indices
    }
}

impl DegreeSetGenerator for RandomTriangularDegreeSet {
    fn degree_set(&mut self, coded_id: usize) -> (Vec<usize>, Vec<usize>) {
        let mut rng = StdRng::seed_from_u64(coded_id as u64);
        let active_indices: Vec<usize>;
        if coded_id < self.k {
            //if self.params.solver_type == SolverType::OrdEnc || self.params.solver_type == SolverType::OrdDec {
            //    panic!("Coded ID {} is out of range for regular LT encoding", coded_id);
            //}
            active_indices = self.triangular_random(&mut rng, coded_id);
        //} //else if coded_id < self.params.num_total() {
        // panic!("Coded ID {} is for parity check constraints", coded_id);
        } else {
            let active_degree = self.active_sampling.random(&mut rng);
            active_indices = sampling_indices_with_degree(active_degree, self.num_active, &mut rng);
        }
        let inactive_degree = self.inactive_sampling.random(&mut rng);
        let inactive_indices =
            sampling_indices_with_degree(inactive_degree, self.num_inactive, &mut rng);
        (active_indices, inactive_indices)
    }
    fn name(&self) -> &str {
        "RandomTriangular"
    }
}
