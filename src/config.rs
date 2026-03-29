// Copyright (c) 2025 Shenghao Yang. All rights reserved.
// Licensed under the MIT License. See LICENSE-MIT for details.

#![allow(deprecated)] // CodeConfig is deprecated; only used by raptor_q to build scheme-compatible configs

use super::degree_sets::DegreeSetConfig;
//use super::precodes::PrecodeConfig;
use fountain_engine::traits::{CodeScheme, HDPC, LDPC};
use fountain_engine::types::{CodeParams, CodeType};

/// Legacy fountain code configuration combining degree sets, precodes, and code parameters.
/// Implemented by `fountain_engine::CodeScheme` for compatibility (e.g. with the `raptor_q` crate).
///
/// **Deprecated:** Use the code scheme types instead: `RandomLTCode`, `LDPCLTCode`, `HDPCLTCode`, `LdpcSysCode`, `HDPCSysCode`.
/// This type is only constructed by `raptor_q::RaptorQConfig` / `RaptorQConfigBuilder::build()`.
#[deprecated(
    since = "1.0.0",
    note = "Use RandomLTCode, LDPCLTCode, HDPCLTCode, LdpcSysCode, or HDPCSysCode instead"
)]
#[derive(Clone, Debug)]
pub struct CodeConfig {
    /// Number of source symbols.
    pub k: usize,
    /// Strategy for generating coded-symbol degree sets.
    pub degree_set_config: DegreeSetConfig,
    /// LDPC and HDPC precode selection.
    pub precode_config: PrecodeConfig,
    /// Ordinary or systematic encoding mode.
    pub code_type: CodeType,
    /// Derived code parameters (k, a, l, h).
    pub params: CodeParams,
}

impl CodeScheme for CodeConfig {
    fn get_params(&self) -> CodeParams {
        self.params.clone()
    }

    fn code_type(&self) -> CodeType {
        self.code_type
    }

    fn create_degree_set_fn(&self) -> Box<dyn FnMut(usize) -> (Vec<usize>, Vec<usize>)> {
        let mut degree_set = self.degree_set_config.create(&self.params);
        Box::new(move |coded_id| degree_set.degree_set(coded_id))
    }

    fn create_precode(&self) -> (Option<Box<dyn HDPC>>, Option<Box<dyn LDPC>>) {
        self.precode_config.generate(&self.params)
    }

    fn max_inactive_num(&self) -> usize {
        (self.params.h + self.params.b + 5).min(128)
    }
}
