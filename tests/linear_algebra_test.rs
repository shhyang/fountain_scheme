// Copyright (c) 2025 Shenghao Yang. All rights reserved.
// Licensed under the MIT License. See LICENSE-MIT for details.

/// Test the linear algebra module together with the lu_solve function in data_manager
/// We us a random m x n matrix A and a sequence of n x t vectors X to test
/// the matrix multiplication B = A X, and the solving of A X = B.
/// The parameters m, n, t can be specified by the user.
/// Then we generate a random m x n matrix A and a sequence of n x t vectors X.
/// We then calculate the matrix multiplication B = A X using both the
/// the linear algebra module and the data manager, and compare the results.
use fountain_engine::algebra::finite_field::GF256;
use fountain_engine::algebra::linear_algebra::Matrix;
use fountain_engine::{DataManager, lu_solve, lu_solve_incr};
use fountain_utility::VecDataOperater;

mod common;
use common::*;

#[test]
fn test_hdpc_matrix_multiplication_31() {
    let field = GF256::default();
    let matrix_d = vec![
        vec![167, 98, 230, 49, 183, 117, 222, 88],
        vec![44, 151, 66, 74, 152, 164, 152, 197],
    ];
    let matrix_x = vec![
        vec![44, 0, 0],
        vec![0, 0, 0],
        vec![0, 0, 0],
        vec![0, 0, 0],
        vec![0, 0, 0],
        vec![1, 0, 0],
        vec![146, 0, 0],
        vec![179, 0, 0],
    ];
    let matrix_b = Matrix::multiply(&field, &matrix_d, &matrix_x);
    for i in 0..matrix_b.len() {
        println!("matrix_b[{}]: {:?}", i, matrix_b[i]);
    }
}

///！ test matrix multiplication
#[test]
fn test_hdpc_matrix_multiplication_21() {
    let field = GF256::default();

    let matrix_d = vec![
        vec![126, 196, 156, 41, 97, 235, 204, 246, 50],
        vec![253, 147, 34, 82, 194, 206, 130, 246, 100],
        vec![225, 60, 69, 165, 159, 134, 29, 246, 201],
        vec![35, 158, 103, 127, 159, 8, 109, 36, 214],
    ];

    let matrix_x = vec![
        vec![0],
        vec![0],
        vec![0],
        vec![0],
        vec![1],
        vec![18],
        vec![36],
        vec![72],
        vec![129],
    ];

    let matrix_b = Matrix::multiply(&field, &matrix_d, &matrix_x);

    for i in 0..matrix_b.len() {
        println!("matrix_b[{}]: {:?}", i, matrix_b[i]);
    }
}

#[test]
fn test_lu_solve_incr() {
    let m = 2;
    let n = 3;
    let t = 3; // vector length
    let field = GF256::default();

    let mut matrix_a = gen_random_matrix(m, n, 1);
    let matrix_x = gen_random_matrix(n, t, 8);
    let matrix_b = Matrix::multiply(&field, &matrix_a, &matrix_x);

    let mut data_manager = DataManager::new_with_operator(Box::new(VecDataOperater::new(t)));
    let mut data_ids = (0..m).collect::<Vec<usize>>();
    for i in 0..m {
        data_manager.insert_data_vector(data_ids[i], &matrix_b[i]);
    }

    let r = 0;
    let mut perm_q = (0..n).collect::<Vec<usize>>();

    let (mut p, mut r) = Matrix::lu_decomp_incr(&field, &mut matrix_a, &mut perm_q, r);

    while r < n {
        matrix_a.truncate(r);
        let mut new_data_ids = Vec::new();
        for i in 0..r {
            new_data_ids.push(data_ids[p[i]]);
        }
        data_ids = new_data_ids;

        let vector_a = gen_random_vector(n, 8);
        let vector_b = Matrix::multiply(&field, &[vector_a.clone()], &matrix_x);
        let data_id = data_manager.temp_data_id();
        data_manager.ensure_zero(&[data_id]);
        data_manager.insert_data_vector(data_id, &vector_b[0]);
        data_ids.push(data_id);
        matrix_a.push(vector_a);
        (p, r) = Matrix::lu_decomp_incr(&field, &mut matrix_a, &mut perm_q, r);
    }

    for i in 0..matrix_a.len() {
        println!("matrix_a[{}]: {:?}", i, matrix_a[i]);
    }

    matrix_a.truncate(r);

    // permute the data ids by the permutation p
    let mut p_data_ids = Vec::new();
    for i in 0..r {
        p_data_ids.push(data_ids[p[i]]);
    }

    // permute the matrix a by the permutation q and run lu_solve
    //let mut matrix_a_new = Vec::new();
    //for i in 0..r {
    //    let mut row = Vec::new();
    //    for j in 0..matrix_a[i].len() {
    //        row.push(matrix_a[i][perm_q[j]]);
    //    }
    //    matrix_a_new.push(row);
    //}
    //lu_solve!(data_manager, &matrix_a_new, &p_data_ids);

    // TODO: Fix lu_solve_incr implementation - currently failing due to permutation handling
    // The basic lu_solve works correctly, but lu_solve_incr needs further investigation
    lu_solve_incr!(data_manager, &matrix_a, &p_data_ids, &perm_q);

    let mut q_data_ids = vec![0; r];
    for i in 0..r {
        q_data_ids[perm_q[i]] = p_data_ids[i];
    }

    let matrix_x_data_manager: Vec<Vec<u8>> = q_data_ids
        .iter()
        .map(|&id| data_manager.get_data_vector(id).to_vec())
        .collect();
    assert_eq!(matrix_x_data_manager, matrix_x);
}

#[test]
fn test_lu_solve() {
    let m = 5;
    let n = 3;
    let t = 4; // vector length
    let field = GF256::default();

    let mut matrix_a = gen_random_matrix(m, n, 8);
    let matrix_x = gen_random_matrix(n, t, 8);
    let matrix_b = Matrix::multiply(&field, &matrix_a, &matrix_x);

    let mut data_manager = DataManager::new_with_operator(Box::new(VecDataOperater::new(t)));
    let data_ids = (0..m).collect::<Vec<usize>>();
    for i in 0..m {
        data_manager.insert_data_vector(data_ids[i], &matrix_b[i]);
    }

    let (p, r) = Matrix::lu_decomp(&field, &mut matrix_a);
    if r != n {
        panic!("The matrix is not invertible, rank = {}", r);
    }

    for i in 0..m {
        println!("matrix_a[{}]: {:?}", i, matrix_a[i]);
    }

    // permute the data ids by the permutation p
    let mut p_data_ids = Vec::new();
    for i in 0..r {
        p_data_ids.push(data_ids[p[i]]);
    }

    matrix_a.truncate(r);

    lu_solve!(data_manager, &matrix_a, &p_data_ids);
    let matrix_x_data_manager: Vec<Vec<u8>> = p_data_ids
        .iter()
        .map(|&id| data_manager.get_data_vector(id).to_vec())
        .collect();
    assert_eq!(matrix_x_data_manager, matrix_x);
}

#[test]
fn test_matrix_multiplication() {
    let m = 10; // number of rows of A
    let n = 4; // number of columns of A
    let t = 3; // vector length
    let field = GF256::default();

    // generate a random m x n matrix A
    let matrix_a = gen_random_matrix(m, n, 8);

    // generate a sequence of n x t vectors X
    let matrix_x = gen_random_matrix(n, t, 8);

    // calculate B = A X using the linear algebra module
    let matrix_b = Matrix::multiply(&field, &matrix_a, &matrix_x);

    let mut data_manager = DataManager::new_with_operator(Box::new(VecDataOperater::new(t)));
    // add the rows of X to the data manager
    for i in 0..n {
        data_manager.insert_data_vector(i, &matrix_x[i]);
    }
    // initialize the result vectors B with zeros
    for i in 0..m {
        data_manager.ensure_zero(&[n + i]);
    }
    // calculate B = A X using the data manager
    for i in 0..m {
        for j in 0..n {
            data_manager.mul_add(j, matrix_a[i][j], n + i);
        }
    }
    let matrix_b_data_manager: Vec<Vec<u8>> = (n..n + m)
        .map(|id| data_manager.get_data_vector(id).to_vec())
        .collect();

    // compare the results
    assert_eq!(matrix_b, matrix_b_data_manager);
}
