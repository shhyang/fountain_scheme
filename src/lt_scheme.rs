// Copyright (c) 2025 Shenghao Yang. All rights reserved.
// Licensed under the MIT License. See LICENSE-MIT for details.

use crate::degree_sets::{
    ideal_soliton_cdf, robust_soliton_cdf, RaptorDegreeSetGenerator, CyclicStrideDegreeSetGenerator,
};
use crate::types::PseudoRandom;

/// *Fountain Code Schemes*
///
/// This module contains the implementation of the fountain code schemes.
use fountain_engine::traits::CodeScheme;
use fountain_engine::types::{CodeParams, CodeType, DecodingConfig};
use crate::types::ProbValue;
use crate::degree_sets::validate_cdf;

/// Random LT Code implementation
///
/// This struct provides a random LT code implementation using ideal soliton
/// degree distribution for efficient fountain coding.
#[derive(Clone)]
pub struct RandomLTCode <R: PseudoRandom> {
    degree_cdf: Vec<ProbValue>,
    params: CodeParams,
    dec_config: DecodingConfig,
    rng: R,
}

impl <R: PseudoRandom> RandomLTCode<R> {
    /// Create a new Random LT Code configuration
    /// # Arguments
    /// * `k` - Number of source symbols
    /// * `degree_distribution` - Degree distribution
    pub fn new(k: usize, cdf: Vec<ProbValue>, rng: R) -> Self {
        let m = cdf.last().expect("CDF must be non-empty");
        if !validate_cdf(&cdf, *m as usize) {
            panic!("Invalid CDF: {:?}", cdf);
        }
        Self {
            degree_cdf: cdf,
            params: CodeParams::new(k, k, 0, 0),
            dec_config: DecodingConfig::default(),
            rng,
        }
    }

    /// LT code with ideal soliton degree distribution (replaces `CodeConfig::lt_ideal`).
    pub fn new_from_ideal_soliton(k: usize, rng: R) -> Self {
        let cdf = ideal_soliton_cdf(k, k);
        Self::new(k, cdf, rng)
    }

    /// LT code with robust soliton degree distribution (replaces `CodeConfig::lt_robust`).
    pub fn new_from_robust_soliton(k: usize, rng: R) -> Self {
        let cdf = robust_soliton_cdf(k, k, 0.1, 0.5);
        Self::new(k, cdf, rng)
    }
}

impl<R: PseudoRandom + 'static> CodeScheme for RandomLTCode<R> {
    /// Get the code parameters
    fn get_params(&self) -> CodeParams {
        self.params.clone()
    }

    /// Get the code type (Ordinary or Systematic)
    fn code_type(&self) -> CodeType {
        CodeType::Ordinary
    }

    /// Create a degree set function for generating coded symbols
    fn create_degree_set_fn(&self) -> Box<dyn FnMut(usize) -> Vec<usize>> {
        let mut degree_set = RaptorDegreeSetGenerator::new(
            self.degree_cdf.clone(),
            self.params.clone(),
            self.rng.clone(),
        );
        Box::new(move |coded_id| degree_set.degree_set(coded_id))
    }

    /// Create precode instances (None for basic LT codes)
    fn create_precode(
        &self,
    ) -> (
        Option<Box<dyn fountain_engine::traits::HDPC>>,
        Option<Box<dyn fountain_engine::traits::LDPC>>,
    ) {
        (None, None)
    }

    /// Get the maximum number of inactive symbols
    fn decoding_config(&self) -> DecodingConfig {
        self.dec_config.clone()
    }
}

/// Deterministic LT Code implementation
///
/// This struct provides a deterministic LT code implementation using a given degree set.
#[derive(Clone)]
pub struct DeterministicLTCode {
    degree_cdf: Vec<ProbValue>,
    params: CodeParams,
}

impl DeterministicLTCode {
    pub fn new(k: usize, cdf: Vec<ProbValue>) -> Self {
        let m = cdf.last().expect("CDF must be non-empty");
        if !validate_cdf(&cdf, *m as usize) {
            panic!("Invalid CDF: {:?}", cdf);
        }
        Self {
            degree_cdf: cdf,
            params: CodeParams::new(k, k, 0, 0),
        }
    }
    /// LT code with ideal soliton degree distribution (replaces `CodeConfig::lt_ideal`).
    pub fn new_from_ideal_soliton(k: usize) -> Self {
        let cdf = ideal_soliton_cdf(k, k);
        Self::new(k, cdf)
    }

    /// LT code with robust soliton degree distribution (replaces `CodeConfig::lt_robust`).
    pub fn new_from_robust_soliton(k: usize) -> Self {
        let cdf = robust_soliton_cdf(k, k, 0.1, 0.5);
        Self::new(k, cdf)
    }
}

impl CodeScheme for DeterministicLTCode {
    fn get_params(&self) -> CodeParams {
        self.params.clone()
    }

    fn code_type(&self) -> CodeType {
        CodeType::Ordinary
    }

    fn create_degree_set_fn(&self) -> Box<dyn FnMut(usize) -> Vec<usize>> {
        let mut degree_set = CyclicStrideDegreeSetGenerator::new(self.params.k, self.degree_cdf.clone());
        Box::new(move |coded_id| degree_set.degree_set(coded_id))
    }

    fn create_precode(
        &self,
    ) -> (
        Option<Box<dyn fountain_engine::traits::HDPC>>,
        Option<Box<dyn fountain_engine::traits::LDPC>>,
    ) {
        (None, None)
    }

    fn decoding_config(&self) -> DecodingConfig {
        DecodingConfig::default()
    }

}

