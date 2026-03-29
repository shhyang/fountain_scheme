// Copyright (c) 2026 Shenghao Yang.
// All rights reserved.

use fountain_engine::traits::HDPC;
use fountain_engine::types::CodeParams;
use fountain_engine::DataManager;
use fountain_utility::VecDataOperater;

/// Validates that `mul_sparse_sh` and `mul_binary` produce consistent results for \(`D_s` `S_h`\).
///
/// This function compares `mul_sparse_sh` and `mul_binary` when calculating \(`D_s` `S_h`\),
/// where \(`S_h`\) is a random binary sparse matrix of size \(l \times h\).
///
/// # Arguments
/// - `hdpc`: An instance of a type implementing the `HDPC` trait
/// - `params`: The coding parameters
///
/// # Panics
/// Panics if the two methods produce different results.
pub fn validate_hdpc_mul_sparse<T: HDPC + ?Sized>(
    hdpc: &T,
    params: &CodeParams,
    sh_fn: &dyn Fn(usize) -> Vec<usize>,
) {
    let h = params.h;
    let l = params.l;
    let a = params.a;

    // Calculate D_s S_h using mul_sparse_sh
    let result_mul_sparse_sh = hdpc.mul_sparse_sh(params, &sh_fn);

    // Build binary column callback v for mul_binary: X has (k+l) rows, same layout as for sparse.
    let v = |row: usize| -> Vec<u8> {
        let mut col = vec![0u8; h];
        if row < a {
            // Rows 0 to a-1: zero (D_a columns)
        } else if row < a + l {
            // Rows a to a+l-1: form S_h (D_s columns)
            for &j in &sh_fn(row - a) {
                col[j] = 1;
            }
        }
        // Rows a+l to a+l+b-1: zero (D_b columns)
        col
    };

    let result_mul_binary = hdpc.mul_binary(params, h, &v);

    // Compare results
    assert_eq!(
        result_mul_sparse_sh.len(),
        h,
        "HDPC validation failed: mul_sparse_sh returned {} rows, expected {}",
        result_mul_sparse_sh.len(),
        h
    );
    assert_eq!(
        result_mul_binary.len(),
        h,
        "HDPC validation failed: mul_binary returned {} rows, expected {}",
        result_mul_binary.len(),
        h
    );

    for i in 0..h {
        assert_eq!(
            result_mul_sparse_sh[i].len(),
            h,
            "HDPC validation failed: mul_sparse_sh row {} has {} columns, expected {}",
            i,
            result_mul_sparse_sh[i].len(),
            h
        );
        assert_eq!(
            result_mul_binary[i].len(),
            h,
            "HDPC validation failed: mul_binary row {} has {} columns, expected {}",
            i,
            result_mul_binary[i].len(),
            h
        );

        for j in 0..h {
            assert_eq!(
                result_mul_sparse_sh[i][j],
                result_mul_binary[i][j],
                "HDPC validation failed: mul_sparse_sh and mul_binary differ at position ({}, {}) for D_s S_h: {} vs {}",
                i,
                j,
                result_mul_sparse_sh[i][j],
                result_mul_binary[i][j]
            );
        }
    }
}

/// Validates that `mul_binary` and `mul_data` produce consistent results.
///
/// This function compares `mul_binary` and `mul_data` for calculating \(DX\)
/// where \(X\) is a binary matrix. This validation requires `fountain_utility`.
///
/// # Arguments
/// - `hdpc`: An instance of a type implementing the `HDPC` trait
/// - `params`: The coding parameters
///
/// # Panics
/// Panics if `mul_binary` and `mul_data` produce different results.
///
#[allow(clippy::too_many_lines)]
pub fn validate_hdpc_mul_data<T: HDPC + ?Sized>(
    hdpc: &T,
    params: &CodeParams,
    matrix_x: &[Vec<u8>],
) {
    let h = params.h;
    let kl = params.num_message_ldpc();
    let vector_length = matrix_x[0].len();

    assert!(kl == matrix_x.len(), "HDPC validation failed: kl != matrix_x.len()");

    // Case I: X is a (k+l)-row matrix
    let v = |row: usize| -> Vec<u8> { matrix_x[row].clone() };

    // Calculate D * X using mul_binary
    let result_mul_binary_dx = hdpc.mul_binary(params, vector_length, &v);

    // Calculate D * X using mul_data
    let mut manager = DataManager::new_with_operator(Box::new(VecDataOperater::new(vector_length)));

    let x_ids: Vec<usize> = (0..kl).collect();
    for i in 0..kl {
        manager.insert_data_vector(x_ids[i], &matrix_x[i]);
    }

    let y_ids: Vec<usize> = (kl..kl + h).collect();
    for &id in &y_ids {
        manager.ensure_zero(&[id]);
    }

    hdpc.mul_data(&mut manager, params, &x_ids, &y_ids);

    let result_mul_data: Vec<Vec<u8>> = y_ids
        .iter()
        .map(|&id| manager.get_data_vector(id).to_vec())
        .collect();

    // Compare results
    for i in 0..h {
        assert_eq!(
            result_mul_binary_dx[i].len(),
            vector_length,
            "HDPC validation failed: mul_binary row {} has {} columns, expected {}",
            i,
            result_mul_binary_dx[i].len(),
            vector_length
        );
        assert_eq!(
            result_mul_data[i].len(),
            vector_length,
            "HDPC validation failed: mul_data row {} has {} columns, expected {}",
            i,
            result_mul_data[i].len(),
            vector_length
        );

        for j in 0..vector_length {
            assert_eq!(
                result_mul_binary_dx[i][j],
                result_mul_data[i][j],
                "HDPC validation failed: mul_binary and mul_data differ at position ({}, {}) for D*X: {} vs {}",
                i,
                j,
                result_mul_binary_dx[i][j],
                result_mul_data[i][j]
            );
        }
    }

    // Case II: X' is the first (a+l)-row matrix of X
    let a = params.a;
    let l = params.l;
    let al = a + l;
    let v_2 = |row: usize| -> Vec<u8> {
        if row < al {
            matrix_x[row].clone()
        } else {
            vec![0u8; vector_length]
        }
    };

    let result_mul_binary_dx_2 = hdpc.mul_binary(params, vector_length, &v_2);

    let mut manager = DataManager::new_with_operator(Box::new(VecDataOperater::new(vector_length)));
    let x_ids: Vec<usize> = (0..al).collect();
    for i in x_ids.clone() {
        manager.insert_data_vector(i, &matrix_x[i]);
    }
    let y_ids: Vec<usize> = (al..al + h).collect();
    for &id in &y_ids {
        manager.ensure_zero(&[id]);
    }
    hdpc.mul_data(&mut manager, params, &x_ids, &y_ids);
    let result_mul_data_2: Vec<Vec<u8>> =
        y_ids.iter().map(|&id| manager.get_data_vector(id).to_vec()).collect();

    // Compare results
    for i in 0..h {
        assert_eq!(
            result_mul_binary_dx_2[i].len(),
            vector_length,
            "HDPC validation failed: mul_binary row {} has {} columns, expected {}",
            i,
            result_mul_binary_dx_2[i].len(),
            vector_length
        );
        assert_eq!(
            result_mul_data_2[i].len(),
            vector_length,
            "HDPC validation failed: mul_data row {} has {} columns, expected {}",
            i,
            result_mul_data_2[i].len(),
            vector_length
        );

        for j in 0..vector_length {
            assert_eq!(
                result_mul_binary_dx_2[i][j],
                result_mul_data_2[i][j],
                "HDPC validation failed: `mul_binary` and `mul_data` differ at position ({}, {}) for D*X': {} vs {}",
                i,
                j,
                result_mul_binary_dx_2[i][j],
                result_mul_data_2[i][j]
            );
        }
    }
}

/// Cross-validates HDPC by running `mul_sparse_sh` vs `mul_binary` and `mul_binary` vs `mul_data` with random inputs.
/// Test-only: compiled only when building tests (`cargo test`).
#[cfg(test)]
pub fn cross_validate_hdpc<T: HDPC + ?Sized>(hdpc: &T, params: &CodeParams) {
    use rand::rngs::StdRng;
    use rand::seq::SliceRandom;
    use rand::{Rng, SeedableRng};

    let h = params.h;
    let kl = params.num_message_ldpc();
    let mut rng = rand::thread_rng();

    let sh_fn = |row: usize| -> Vec<usize> {
        let mut r = StdRng::seed_from_u64(row as u64);
        let mut cols: Vec<usize> = (0..h).collect();
        cols.shuffle(&mut r);
        cols.truncate(3.min(h));
        cols
    };
    validate_hdpc_mul_sparse(hdpc, params, &sh_fn);

    let vector_length = 4.max(h / 2);
    let mut matrix_x = vec![vec![0u8; vector_length]; kl];
    for i in 0..kl {
        let num_ones = rng.gen_range(1..=vector_length);
        let mut cols: Vec<usize> = (0..vector_length).collect();
        cols.shuffle(&mut rng);
        for j in 0..num_ones {
            matrix_x[i][cols[j]] = 1;
        }
    }
    validate_hdpc_mul_data(hdpc, params, &matrix_x);
}

