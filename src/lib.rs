// Copyright (c) 2025 Shenghao Yang. All rights reserved.
// Licensed under the MIT License. See LICENSE-MIT for details.

//! Fountain Code Implementation Library
//!
//! This library provides specific implementations of fountain codes and configurations.
//!
//! Prefer the code scheme types (`RandomLTCode`, `LDPCLTCode`, `HDPCLTCode`, `BinaryHDPCLTCode`).
//! Use [`CodeType::Systematic`](fountain_engine::types::CodeType::Systematic) via `as_systematic()` where provided.

/// Degree-set generation strategies and distributions.
pub mod degree_sets;
/// Encoder/decoder for LT codes with HDPC precode.
pub mod hdpc_lt_scheme;
/// Systematic encoder/decoder for LT codes with HDPC precode.
pub mod binary_hdpc_lt_scheme;
/// Encoder/decoder for LT codes with LDPC precode.
pub mod ldpc_lt_scheme;
/// Pure LT code encoder/decoder (no precodes).
pub mod lt_scheme;
/// Parameter-generation strategies that derive `CodeParams` from `k`.
pub mod parameters;
/// LDPC and HDPC precode implementations.
pub mod precodes;
/// Internal integer math helpers shared by modules.
pub(crate) mod math_utils;

pub mod types;

/// Cross-validation utilities for LDPC and HDPC implementations.
pub mod validation;

pub use hdpc_lt_scheme::*;
pub use binary_hdpc_lt_scheme::*;
pub use ldpc_lt_scheme::*;
pub use lt_scheme::*;

pub use validation::*;
