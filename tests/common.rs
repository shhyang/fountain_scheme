// Copyright (c) 2025 Shenghao Yang. All rights reserved.
// Licensed under the MIT License. See LICENSE-MIT for details.

use fountain_engine::DataManager;
use fountain_engine::algebra::finite_field::GF256;
use fountain_utility::VecDataOperater;
use rand::Rng;

/// In-memory [`DataManager`] whose operator uses the same GF(256) as `field`.
#[allow(dead_code)]
pub fn data_manager_for_gf256(vector_len: usize, field: &GF256) -> DataManager {
    let mut manager = DataManager::new_with_operator(Box::new(VecDataOperater::new_with_gf256(
        vector_len,
        field.primitive_polynomial(),
    )));
    manager.config_finite_field_from(field);
    manager
}
/// Generate a random m x n matrix over the finite field GF(2^ffsize)
#[must_use]
pub fn gen_random_matrix(m: usize, n: usize, ffsize: u8) -> Vec<Vec<u8>> {
    let mut matrix_a = Vec::new();
    // All arithmetic in `u32` so we never cast from signed `i32` to `u8`.
    let shift = u32::from(ffsize.saturating_sub(1));
    let q = u8::try_from((2u32 << shift).saturating_sub(1)).expect("ffsize too large for u8");
    for _ in 0..m {
        let mut row = Vec::new();
        for _ in 0..n {
            row.push(rand::thread_rng().gen_range(0..=q));
        }
        matrix_a.push(row);
    }
    matrix_a
}

/// Generate a random length-`n` vector over GF(2^`ffsize`).
#[allow(dead_code)] // only some integration test crates use this helper
#[must_use]
pub fn gen_random_vector(n: usize, ffsize: u8) -> Vec<u8> {
    let shift = u32::from(ffsize.saturating_sub(1));
    let q = u8::try_from((2u32 << shift).saturating_sub(1)).expect("ffsize too large for u8");
    (0..n)
        .map(|_| rand::thread_rng().gen_range(0..=q))
        .collect()
}
