// Copyright (c) 2026 Shenghao Yang.
// All rights reserved.

//! Cross-validation utilities for verifying LDPC and HDPC precode consistency.

/// HDPC `mul_binary` / `mul_sparse_sh` / `mul_data` consistency checks.
pub mod hdpc;
/// LDPC row-vs-column consistency checks.
pub mod ldpc;
pub mod pseudo_rand;

pub use hdpc::*;
pub use ldpc::*;

//pub use rand::*;
