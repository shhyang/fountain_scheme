// Copyright (c) 2025 Shenghao Yang. All rights reserved.
// Licensed under the MIT License. See LICENSE-MIT for details.

use fountain_engine::types::CodeParams;

/// Predefined parameter-generation strategies that compute `CodeParams` from `k`.
#[derive(Clone, Debug)]
pub enum ParamsConfig {
    /// Pure LT code: `l = h = 0`, all symbols are active.
    LT,
    /// LT code with LDPC precode: sets `l` based on `k`.
    LDPCLT,
    /// LDPC-LT with pre-inactivation: sets both `l` and reduces `a` by an inactivity count.
    LDPCLTPreInactive,
    /// Full Raptor code with LDPC + HDPC: sets `l`, `h`, and pre-inactivation.
    HDPC,
    /// Lookup a parameter strategy by uppercase name at runtime.
    ByName(String),
}

impl ParamsConfig {
    /// Computes the [`CodeParams`] for `k` source symbols using this strategy.
    pub fn generate(&self, k: usize) -> CodeParams {
        match self {
            ParamsConfig::LT => LTParams.generate(k),
            ParamsConfig::LDPCLT => LDPCLTParams.generate(k),
            ParamsConfig::LDPCLTPreInactive => LDPCLTParamsPreInactive.generate(k),
            ParamsConfig::HDPC => HDPCParams.generate(k),
            ParamsConfig::ByName(name) => create_params_by_name(name, k),
        }
    }
}

/// Create parameters by name using the ParamsGenerator trait
pub fn create_params_by_name(name: &str, k: usize) -> CodeParams {
    match name.to_uppercase().as_str() {
        "LT" => LTParams.generate(k),
        "LDPC_LT" => LDPCLTParams.generate(k),
        "LDPC_LT_PRE_INACTIVE" => LDPCLTParamsPreInactive.generate(k),
        "HDPC" => HDPCParams.generate(k),
        //"CUSTOM" => CustomParams.generate(k),
        _ => panic!("Unknown parameter generator: {}", name),
    }
}

/*
#[derive(Clone, Debug)]
pub enum ParamsGeneratorEnum {
    LT, // LT code, no precode
    LDPC_LT, // LT code with LDPC precode
    LDPC_LT_PRE_INACTIVE, // LT code with LDPC precode and pre-inactivation.
    HDPC, // LDPC and HDPC precode
    //Custom(String),
}

impl ParamsGeneratorEnum {
    pub fn generate(&self, k: usize) -> CodeParams {
        match self {
            ParamsGeneratorEnum::LT => LTParams.generate(k),
            ParamsGeneratorEnum::LDPC_LT => LDPCLTParams.generate(k),
            ParamsGeneratorEnum::LDPC_LT_PRE_INACTIVE => LDPCLTParamsPreInactive.generate(k),
            ParamsGeneratorEnum::HDPC => HDPCParams.generate(k),
        }
    }
}
*/

/// Trait for types that derive [`CodeParams`] from a source-symbol count `k`.
pub trait ParamsGenerator {
    /// Compute the code parameters for `k` source symbols.
    fn generate(&self, k: usize) -> CodeParams;
    /// Human-readable name for this parameter strategy.
    fn name(&self) -> &str;
}

/// Pure LT parameters: no LDPC, no HDPC, all symbols active.
pub struct LTParams;
impl ParamsGenerator for LTParams {
    fn generate(&self, k: usize) -> CodeParams {
        CodeParams::new(k, k, 0, 0)
    }
    fn name(&self) -> &str {
        "LT"
    }
}

/// LT parameters with LDPC precode rows, `h = 0`.
pub struct LDPCLTParams;
impl LDPCLTParams {
    fn ldpc_num(k: usize) -> usize {
        match k {
            4..=10 => 1,
            11..=20 => 2,
            21..=30 => 3,
            31..=50 => 5,
            51..=70 => 7,
            71..=110 => 11,
            111..=130 => 13,
            131..=170 => 17,
            _ => (0.1 * k as f64).floor() as usize,
        }
    }
}

impl ParamsGenerator for LDPCLTParams {
    fn generate(&self, k: usize) -> CodeParams {
        CodeParams::new(k, k, Self::ldpc_num(k), 0)
    }
    fn name(&self) -> &str {
        "LDPC_LT"
    }
}

/// LDPC-LT parameters with pre-inactivation (reduces active count `a`).
pub struct LDPCLTParamsPreInactive;

impl LDPCLTParamsPreInactive {
    fn inactive_num(k: usize) -> usize {
        match k {
            0..=3 => 1,
            4..=7 => 1,
            8..=15 => 2,
            16..=31 => 4,
            32..=63 => 5,
            64..=127 => 6,
            _ => (k as f64).sqrt() as usize,
        }
    }

    fn ldpc_num(k: usize) -> usize {
        match k {
            4..=20 => 2,
            21..=30 => 3,
            31..=50 => 5,
            51..=70 => 7,
            71..=110 => 11,
            111..=130 => 13,
            131..=170 => 17,
            _ => (0.1 * k as f64).floor() as usize,
        }
    }
}

impl ParamsGenerator for LDPCLTParamsPreInactive {
    fn generate(&self, k: usize) -> CodeParams {
        CodeParams::new(k, k - Self::inactive_num(k), Self::ldpc_num(k), 0)
    }
    fn name(&self) -> &str {
        "LDPC_LT_PRE_INACTIVE"
    }
}

/// Full Raptor parameters with LDPC rows, HDPC rows, and pre-inactivation.
pub struct HDPCParams;

impl HDPCParams {
    fn inactive_num(k: usize) -> usize {
        match k {
            0..=3 => 1,
            4..=7 => 1,
            8..=15 => 2,
            16..=31 => 4,
            32..=63 => 5,
            64..=127 => 6,
            _ => (k as f64).sqrt() as usize,
        }
    }

    fn ldpc_num(k: usize) -> usize {
        match k {
            4..=10 => 1,
            11..=20 => 2,
            21..=30 => 3,
            31..=50 => 5,
            51..=70 => 7,
            71..=110 => 11,
            111..=130 => 13,
            131..=170 => 17,
            _ => (0.1 * k as f64).floor() as usize,
        }
    }

    fn hdpc_num(k: usize) -> usize {
        match k {
            1..=20 => 2,
            21..=30 => 4,
            31..=50 => 6,
            51..=70 => 8,
            _ => 10,
        }
    }
}

impl ParamsGenerator for HDPCParams {
    fn generate(&self, k: usize) -> CodeParams {
        CodeParams::new(
            k,
            k - Self::inactive_num(k),
            Self::ldpc_num(k),
            Self::hdpc_num(k),
        )
    }
    fn name(&self) -> &str {
        "HDPC"
    }
}
