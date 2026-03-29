// Copyright (c) 2026 Shenghao Yang.
// All rights reserved.

//! Cross-validation utilities for verifying LDPC and HDPC precode consistency.

/// LDPC row-vs-column consistency checks.
pub mod ldpc;
/// HDPC `mul_binary` / `mul_sparse_sh` / `mul_data` consistency checks.
pub mod hdpc;
pub mod pseudo_rand;

pub use ldpc::*;
pub use hdpc::*;

//pub use rand::*;