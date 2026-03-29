// Copyright (c) 2025 Shenghao Yang. All rights reserved.
// Licensed under the MIT License. See LICENSE-MIT for details.

use super::{default_HDPC, NoHDPC, NoLDPC, RQLDPC, ReversedLDPC, R10LDPC, R10HDPC};
use fountain_engine::CodeParams;
use fountain_engine::traits::HDPC;
use fountain_engine::traits::LDPC;

/*
#[derive(Clone, Debug)]
pub enum PrecodeType {
    None, // No precode
    LDPC, // RQ LDPC only precode
    ReversedLDPC, // Reversed LDPC only precode
    LDPC_HDPC, // RQ LDPC and HDPC precode
    ReversedLDPC_HDPC, // Reversed LDPC and RQHDPC precode, used for systematic encoding
}

impl PrecodeType {
    pub fn generate_precode(&self, params: CodeParams) -> (Option<DefaultHDPC>, LDPCType) {
        let (hdpc, ldpc) = match self {
            PrecodeType::None => (None, LDPCType::None),
            PrecodeType::LDPC => {
                if params.l == 0 {
                    (None, LDPCType::None)
                } else {
                    (None, LDPCType::RQLDPC)
                }
            },
            PrecodeType::ReversedLDPC => {
                if params.l == 0 {
                    (None, LDPCType::None)
                } else {
                    (None, LDPCType::ReversedLDPC)
                }
            },
            PrecodeType::LDPC_HDPC => {
                let ldpc = if params.l == 0 {
                    LDPCType::None
                } else {
                    LDPCType::RQLDPC
                };
                let hdpc = if params.h == 0 {
                    None
                } else {
                    Some(DefaultHDPC::new(params.h))
                };
                (hdpc, ldpc)
            },
            PrecodeType::ReversedLDPC_HDPC => {
                let ldpc = if params.l == 0 {
                    LDPCType::None
                } else {
                    LDPCType::ReversedLDPC
                };
                let hdpc = if params.h == 0 {
                    None
                } else {
                    Some(DefaultHDPC::new(params.h))
                };
                (hdpc, ldpc)
            },
        };
        (hdpc, ldpc)
    }
}

*/

/// Selects which LDPC precode variant to use.
#[derive(Clone, Debug)]
pub enum LDPCType {
    /// No LDPC precode.
    NoLDPC,
    /// Raptor-10 LDPC precode.
    R10LDPC,
    /// Standard RaptorQ LDPC precode (forward column order).
    RQLDPC,
    /// LDPC precode with reversed column order, used for systematic encoding.
    ReversedLDPC,
}

impl LDPCType {
    /// Instantiates the selected LDPC variant from the given code parameters.
    pub fn create(&self, params: &CodeParams) -> Box<dyn LDPC> {
        match self {
            LDPCType::NoLDPC => Box::new(NoLDPC),
            LDPCType::R10LDPC => Box::new(R10LDPC::new(params)),
            LDPCType::RQLDPC => Box::new(RQLDPC::new(params)),
            LDPCType::ReversedLDPC => Box::new(ReversedLDPC::new(params)),
        }
    }
}

/// Selects which HDPC precode variant to use.
#[derive(Clone, Debug)]
pub enum HDPCType {
    /// No HDPC precode.
    NoHDPC,
    /// Raptor-10 HDPC precode.
    R10HDPC,
    /// Deterministic HDPC with cyclic delta column mapping.
    Default,
}

impl HDPCType {
    /// Create an HDPC instance from the type
    pub fn create(&self, params: &CodeParams) -> Box<dyn HDPC> {
        match self {
            HDPCType::NoHDPC => Box::new(NoHDPC),
            HDPCType::R10HDPC => Box::new(R10HDPC::new()),
            HDPCType::Default => default_HDPC(params.h),
        }
    }
}


/// Combined LDPC and HDPC precode configuration.
#[deprecated(since = "1.0.0", note = "Use HDPCType and LDPCType instead")]
#[derive(Clone, Debug)]
pub struct PrecodeConfig {
    /// Which HDPC variant to use (or none).
    pub hdpc_type: HDPCType,
    /// Which LDPC variant to use (or none).
    pub ldpc_type: LDPCType,
}

impl PrecodeConfig {
    /// Constructs boxed HDPC and LDPC trait objects according to this configuration.
    ///
    /// Returns `None` for either precode when the corresponding parameter count is zero.
    pub fn generate(&self, params: &CodeParams) -> (Option<Box<dyn HDPC>>, Option<Box<dyn LDPC>>) {
        let ldpc_option = if params.l == 0 {
            None
        } else {
            Some(self.ldpc_type.create(params))
        };

        let hdpc_option = if params.h == 0 {
            None
        } else {
            Some(self.hdpc_type.create(params))
        };

        (hdpc_option, ldpc_option)
    }
}


