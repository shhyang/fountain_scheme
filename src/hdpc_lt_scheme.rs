// Copyright (c) 2025 Shenghao Yang. All rights reserved.
// Licensed under the MIT License. See LICENSE-MIT for details.

use crate::degree_sets::asymp_optimal_cdf;
use crate::degree_sets::validate_cdf;
use crate::degree_sets::RaptorDegreeSetGenerator;
use crate::parameters::rq_like;
use crate::precodes::default_rq_hdpc;
use crate::precodes::RQLDPC;
use crate::types::ProbValue;
use crate::types::PseudoRandom;

/// *HDPC-LT Code Implementation*
///
/// LT fountain coding with R10-style LDPC + HDPC precoding (full R10 intermediate structure),
/// Raptor-style degree sets, and [`crate::parameters::r10`] parameters.
use fountain_engine::traits::{CodeScheme, HDPC, LDPC};
use fountain_engine::types::{CodeParams, CodeType, DecodingConfig, DegreeSetFn};

/// HDPC-LT code: R10 parameters (nonzero `l` and `h`), [`R10LDPC`], and default cyclic HDPC.
#[derive(Clone)]
pub struct HDPCLTCode<R: PseudoRandom> {
    degree_cdf: Vec<ProbValue>,
    params: CodeParams,
    rng: R,
    code_type: CodeType,
}

fn hdpc_lt_code_params(k: usize) -> CodeParams {
    rq_like::generate_rq_like_parameters(k)
}

impl<R: PseudoRandom> HDPCLTCode<R> {
    /// New HDPC-LT configuration with an explicit degree CDF.
    pub fn new(k: usize, cdf: Vec<ProbValue>, rng: R) -> Self {
        let m = cdf.last().expect("CDF must be non-empty");
        if !validate_cdf(&cdf, *m as usize) {
            panic!("Invalid CDF: {:?}", cdf);
        }
        Self {
            degree_cdf: cdf,
            params: hdpc_lt_code_params(k),
            rng,
            code_type: CodeType::Ordinary,
        }
    }

    /// Convenience: asymptotically optimal CDF (same helper family as [`LDPCLTCode::new_with_ideal_soliton`](crate::LDPCLTCode::new_with_ideal_soliton)).
    pub fn new_with_ideal_soliton(k: usize, rng: R) -> Self {
        let params = hdpc_lt_code_params(k);
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

impl<R: PseudoRandom + 'static> CodeScheme for HDPCLTCode<R> {
    fn get_params(&self) -> CodeParams {
        self.params.clone()
    }

    fn code_type(&self) -> CodeType {
        self.code_type
    }

    fn create_degree_set_fn(&self) -> DegreeSetFn {
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

    fn create_precode(&self) -> (Option<Box<dyn HDPC>>, Option<Box<dyn LDPC>>) {
        let hdpc = (self.params.h > 0).then(|| default_rq_hdpc(self.params.h));
        let ldpc = (self.params.l > 0).then(|| Box::new(RQLDPC::new(&self.params)) as Box<dyn LDPC>);
        (hdpc, ldpc)
    }

    fn decoding_config(&self) -> DecodingConfig {
        DecodingConfig::default().with_max_inact_num((self.params.h + self.params.b + 5).min(128))
    }
}
