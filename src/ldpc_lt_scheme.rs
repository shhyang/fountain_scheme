// Copyright (c) 2025 Shenghao Yang. All rights reserved.
// Licensed under the MIT License. See LICENSE-MIT for details.

use crate::types::ProbValue;
use crate::types::PseudoRandom;
use crate::degree_sets::validate_cdf;
use crate::degree_sets::RaptorDegreeSetGenerator;
use crate::precodes::R10LDPC;
use crate::parameters::r10;
use crate::degree_sets::asymp_optimal_cdf;

/// *LDPC-LT Code Implementation*
///
/// This module contains the implementation of LT codes with LDPC precoding.
/// LDPC precoding provides structured redundancy that improves the performance
/// of fountain codes by adding low-density parity check constraints.
use fountain_engine::traits::{CodeScheme, HDPC, LDPC};
use fountain_engine::types::{CodeParams, CodeType, DecodingConfig, DegreeSetFn};

/// LDPC-LT Code implementation
///
/// This struct provides an LT code implementation with 
/// - R10-LDPC precoding,
/// - Raptor degree set generator, and 
/// - R10 parameters.
#[derive(Clone)]
pub struct LDPCLTCode <R: PseudoRandom> {
    degree_cdf: Vec<ProbValue>,
    params: CodeParams,
    rng: R,
    code_type: CodeType,
}

/// R10-style `k`, `a`, `l` from [`r10::generate_r10_parameters`], with `h = 0`.
///
/// This type only installs LDPC in [`CodeScheme::create_precode`]; leaving `h > 0` would reserve
/// HDPC variables that [`fountain_engine::core::precode::precode_encode`] never fills.
fn ldpc_lt_code_params(k: usize) -> CodeParams {
    let mut p = r10::generate_r10_parameters(k);
    p.h = 0;
    p
}

impl <R: PseudoRandom> LDPCLTCode<R> {
    /// Create a new LDPC-LT Code configuration
    ///
    /// # Arguments
    ///
    /// * `k` - Number of source symbols
    /// * `cdf` - Degree distribution CDF
    /// * `rng` - Random number generator
    ///
    /// # Returns
    ///
    /// A new `LDPCLTCode` instance with LDPC precoding
    pub fn new(k: usize, cdf: Vec<ProbValue>, rng: R) -> Self {
        let m = cdf.last().expect("CDF must be non-empty");
        if !validate_cdf(&cdf, *m as usize) {
            panic!("Invalid CDF: {:?}", cdf);
        }
        Self {
            degree_cdf: cdf,
            params: ldpc_lt_code_params(k),
            rng,
            code_type: CodeType::Ordinary,
        }
    }

    pub fn new_with_ideal_soliton(k: usize, rng: R) -> Self {
        let params = ldpc_lt_code_params(k);
        let cdf = asymp_optimal_cdf(k, 2, k);
        Self {
            degree_cdf: cdf,
            params,
            rng,
            code_type: CodeType::Ordinary,
        }        
    }

    pub fn as_systematic(&mut self) {
        self.code_type = CodeType::Systematic;
    }
}

impl <R: PseudoRandom + 'static> CodeScheme for LDPCLTCode<R> {
    /// Get the code parameters
    fn get_params(&self) -> CodeParams {
        self.params.clone()
    }

    /// Get the code type (Ordinary or Systematic)
    fn code_type(&self) -> CodeType {
        self.code_type
    }

    /// Create a degree set function for generating coded symbols
    fn create_degree_set_fn(&self) -> DegreeSetFn {
        // LT edges attach to message + LDPC variables (`0..num_active`); HDPC indices are separate.
        let mut degree_set = if self.code_type == CodeType::Systematic {
            RaptorDegreeSetGenerator::new_systematic(
                self.degree_cdf.clone(),
                self.params.clone(),
                self.rng.clone(),
            )
        } else {
            RaptorDegreeSetGenerator::new(
                self.degree_cdf.clone(),
                self.params.clone(),
                self.rng.clone(),
            )
        };
        Box::new(move |coded_id| degree_set.degree_set(coded_id))
    }

    /// Create precode instances (LDPC for LDPC-LT codes)
    fn create_precode(&self) -> (Option<Box<dyn HDPC>>, Option<Box<dyn LDPC>>) {
        (None, Some(Box::new(R10LDPC::new(&self.params))))
    }

    /// Get the decoding configuration
    fn decoding_config(&self) -> DecodingConfig {
        DecodingConfig::default().with_max_inact_num((self.params.h + self.params.b + 5).min(128))
    }
}
