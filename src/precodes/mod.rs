// Copyright (c) 2025 Shenghao Yang. All rights reserved.
// Licensed under the MIT License. See LICENSE-MIT for details.

//! LDPC and HDPC precode implementations.

mod gray;
pub mod hdpc;
pub mod ldpc;
pub mod hdpc_binary;

pub use hdpc::*;
pub use ldpc::*;
pub use hdpc_binary::*;
