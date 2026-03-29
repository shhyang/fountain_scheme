// Copyright (c) 2025 Shenghao Yang. All rights reserved.
// Licensed under the MIT License. See LICENSE-MIT for details.

use super::SolitonDistribution;
use super::sampling::*;
use super::traits::DegreeSetGenerator;
use fountain_engine::CodeParams;
use rand::{SeedableRng, rngs::StdRng};

/// Degree set using random sampling.
pub struct RandomDegreeSet {
    active_sampling: DegreeSampling,
    inactive_sampling: DegreeSampling,
    //params: CodeParams, // range of indices for sampling
    num_active: usize,
    num_inactive: usize,
}

impl RandomDegreeSet {
    /// Creates a new random degree set from the given distribution and code parameters.
    pub fn new(degree_distribution: &SolitonDistribution, params: &CodeParams) -> Self {
        let mut active_pdf = degree_distribution.pdf();
        shorten_pdf(&mut active_pdf, params.num_active());
        let active_cdf = pdf_to_cdf(&active_pdf);

        //dbg!("active_cdf", active_pdf.len(), &active_cdf);
        // the degree of inactive packets is 2 or 3
        let mut inactive_pdf = vec![0.0, 0.0, 0.5, 0.5];
        shorten_pdf(&mut inactive_pdf, params.num_inactive());
        let inactive_cdf = pdf_to_cdf(&inactive_pdf);

        //dbg!("inactive_cdf", inactive_pdf.len(), &inactive_cdf);

        let active_sampling = DegreeSampling::new(active_cdf);
        let inactive_sampling = DegreeSampling::new(inactive_cdf);
        Self {
            active_sampling,
            inactive_sampling,
            num_active: params.num_active(),
            num_inactive: params.num_inactive(),
        }
    }
}

impl DegreeSetGenerator for RandomDegreeSet {
    fn degree_set(&mut self, coded_id: usize) -> (Vec<usize>, Vec<usize>) {
        let mut rng = StdRng::seed_from_u64(coded_id as u64);
        let active_degree = self.active_sampling.random(&mut rng);
        let active_indices = sampling_indices_with_degree(active_degree, self.num_active, &mut rng);
        let inactive_degree = self.inactive_sampling.random(&mut rng);
        let inactive_indices =
            sampling_indices_with_degree(inactive_degree, self.num_inactive, &mut rng);
        (active_indices, inactive_indices)
    }
    fn name(&self) -> &str {
        "Random"
    }
}
