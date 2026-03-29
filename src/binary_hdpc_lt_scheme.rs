// Copyright (c) 2025 Shenghao Yang. All rights reserved.
// Licensed under the MIT License. See LICENSE-MIT for details.

use crate::degree_sets::asymp_optimal_cdf;
use crate::degree_sets::validate_cdf;
use crate::degree_sets::RaptorDegreeSetGenerator;
use crate::parameters::rq_like;
use crate::precodes::R10HDPC;
use crate::precodes::R10LDPC;
use crate::types::ProbValue;
use crate::types::PseudoRandom;

/// *Binary HDPC-LT code (R10 LDPC + R10 HDPC)*
///
/// Same R10-style parameters and degree-set construction as [`HDPCLTCode`](crate::HDPCLTCode),
/// but uses [`R10LDPC`](crate::precodes::R10LDPC) / [`R10HDPC`](crate::precodes::R10HDPC).
/// Defaults to [`CodeType::Ordinary`] (same as [`HDPCLTCode`](crate::HDPCLTCode)); call [`Self::as_systematic`]
/// for systematic degree sets and systematic encoder/decoder behavior.
use fountain_engine::traits::{CodeScheme, HDPC, LDPC};
use fountain_engine::types::{CodeParams, CodeType, DecodingConfig, DegreeSetFn};

/// R10 `CodeParams`, [`R10LDPC`], [`R10HDPC`], Raptor-style degree sets (ordinary or systematic).
#[derive(Clone)]
pub struct BinaryHDPCLTCode<R: PseudoRandom> {
    degree_cdf: Vec<ProbValue>,
    params: CodeParams,
    rng: R,
    code_type: CodeType,
}

fn hdpc_sys_code_params(k: usize) -> CodeParams {
    rq_like::generate_rq_like_parameters(k)
}

impl<R: PseudoRandom> BinaryHDPCLTCode<R> {
    /// New configuration with an explicit degree CDF ([`CodeType::Ordinary`] until [`Self::as_systematic`]).
    pub fn new(k: usize, cdf: Vec<ProbValue>, rng: R) -> Self {
        let m = cdf.last().expect("CDF must be non-empty");
        if !validate_cdf(&cdf, *m as usize) {
            panic!("Invalid CDF: {:?}", cdf);
        }
        Self {
            degree_cdf: cdf,
            params: hdpc_sys_code_params(k),
            rng,
            code_type: CodeType::Ordinary,
        }
    }

    /// Convenience: same CDF construction as [`HDPCLTCode::new_with_ideal_soliton`](crate::HDPCLTCode::new_with_ideal_soliton).
    pub fn new_with_ideal_soliton(k: usize, rng: R) -> Self {
        let params = hdpc_sys_code_params(k);
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

impl<R: PseudoRandom + 'static> CodeScheme for BinaryHDPCLTCode<R> {
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
        let hdpc = (self.params.h > 0).then(|| Box::new(R10HDPC::new()) as Box<dyn HDPC>);
        let ldpc = (self.params.l > 0).then(|| Box::new(R10LDPC::new(&self.params)) as Box<dyn LDPC>);
        (hdpc, ldpc)
    }

    fn decoding_config(&self) -> DecodingConfig {
        DecodingConfig::default().with_max_inact_num((self.params.h + self.params.b + 5).min(128))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::validation::pseudo_rand::XorShift64;
    use fountain_engine::types::CodeType;

    #[test]
    fn test_hdpc_sys_creation() {
        let k = 50;
        let cdf = asymp_optimal_cdf(k, 2, k);
        let code = BinaryHDPCLTCode::new(k, cdf, XorShift64::new(1));

        assert_eq!(code.get_params().k, k);
        assert!(code.get_params().h > 0);
        assert!(code.get_params().l > 0);
        assert_eq!(code.code_type(), CodeType::Ordinary);
    }

    #[test]
    fn test_hdpc_sys_with_default_setting() {
        let k = 100;
        let code = BinaryHDPCLTCode::new_with_ideal_soliton(k, XorShift64::new(2));

        assert_eq!(code.get_params().k, k);
        assert!(code.get_params().h > 0);
        assert!(code.get_params().l > 0);
        assert_eq!(code.code_type(), CodeType::Ordinary);
    }

    #[test]
    fn test_degree_set_function() {
        let k = 20;
        let cdf = asymp_optimal_cdf(k, 2, k);
        let code = BinaryHDPCLTCode::new(k, cdf, XorShift64::new(3));
        let mut degree_set_fn = code.create_degree_set_fn();
        let params = code.get_params();
        let set = degree_set_fn(params.num_total());
        assert!(!set.is_empty());
        assert!(set.len() <= params.num_active());
    }
}
