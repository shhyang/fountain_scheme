// Copyright (c) 2025 Shenghao Yang. All rights reserved.
// Licensed under the MIT License. See LICENSE-MIT for details.

use crate::precodes::gray::Gray;
use fountain_engine::DataManager;
use fountain_engine::traits::HDPC;
use fountain_engine::types::CodeParams;
use fountain_engine::algebra::binary_vector::BinaryVector;

fn r10_mul_data(
    manager: &mut DataManager,
    h: usize,
    x_ids: &[usize],
    y_ids: &[usize],
) {
    let m = x_ids.len();
    let weight = if h % 2 == 0 { h / 2 } else { h.div_ceil(2) };
    let mut gray = Gray::new(h, weight);
    gray.with_column(m-1);
    //let kl = params.num_message_ldpc();
    let temp_id = manager.temp_data_id();
    //dbg!(&m, &h);
    //manager.ensure_zero(&[temp_id]);
    manager.add_to_vector(&[x_ids[m-1]], temp_id);
    for i in (0..m - 1).rev() {
        //dbg!("mul_vector", &temp_id, manager.get_vector(temp_id));
        let j = gray.previous_delta()
            .iter()
            .map(|&id| y_ids[id])
            .collect::<Vec<_>>();
        manager.broadcast_add(temp_id, &j);
        manager.add_to_vector(&[x_ids[i]], temp_id);
    }
    //dbg!("mul_vector", &temp_id, manager.get_vector(temp_id));
    //let j = (0..weight).map(|i| y_ids[i]).collect::<Vec<_>>();
    let j = gray.current_column_positions().iter().map(|&i| y_ids[i]).collect::<Vec<_>>();
    manager.broadcast_add(temp_id, &j);
    manager.remove(temp_id);
}

/// Raptor-10 HDPC precode using balanced Gray code binary matrices.
///
/// Multiplication is performed over GF(2) using constant-weight Gray code
/// columns, which allow efficient incremental updates via two-bit deltas.
pub struct R10HDPC;

impl Default for R10HDPC {
    fn default() -> Self {
        Self::new()
    }
}

impl R10HDPC {
    /// Creates a new Raptor-10 binary HDPC precode instance.
    #[must_use]
    pub fn new() -> Self {
        Self
    }
}

impl HDPC for R10HDPC {
    fn mul_data(
        &self,
        manager: &mut DataManager,
        params: &CodeParams,
        x_ids: &[usize],
        y_ids: &[usize],
    ) {
        r10_mul_data(manager, params.h, x_ids, y_ids);
    }

    fn mul_binary(
        &self,
        params: &CodeParams,
        n: usize,
        v: &dyn Fn(usize) -> Vec<u8>,
    ) -> Vec<Vec<u8>> {
        let h = params.h;
        let kl = params.num_message_ldpc();
        let weight = if h % 2 == 0 { h / 2 } else { h.div_ceil(2) };
        let mut gray = Gray::new(h, weight);
        gray.with_column(kl - 1);

        let mut result = vec![vec![0u8; n]; h];

        let mut tmp = vec![0u8; n];
        for (col, &b) in v(kl - 1).iter().enumerate() {
            if b == 1 {
                tmp[col] ^= 1;
            }
        }

        for i in (0..kl - 1).rev() {
            let j = gray.previous_delta();
            BinaryVector::add_inplace(&mut result[j[0]], &tmp);
            BinaryVector::add_inplace(&mut result[j[1]], &tmp);
            for (col, &b) in v(i).iter().enumerate() {
                if b == 1 {
                    tmp[col] ^= 1;
                }
            }
        }
        let j = gray.current_column_positions();
        for i in j {
            BinaryVector::add_inplace(&mut result[i], &tmp);
        }
        result
    }
 
    fn mul_sparse(
        &self,
        params: &CodeParams,
        n: usize,
        s: &dyn Fn(usize) -> Vec<usize>,
    ) -> Vec<Vec<u8>> {
        let v = |i: usize| {
            let mut col = vec![0u8; n];
            for &j in &s(i) {
                col[j] = 1;
            }
            col
        };
        self.mul_binary(params, n, &v)
    }

    fn mul_sparse_sh(
        &self,
        params: &CodeParams,
        s: &dyn Fn(usize) -> Vec<usize>,
    ) -> Vec<Vec<u8>> {
        let h = params.h;
        let l = params.l;
        let a = params.a;
        let al = a + l;
        let weight = if h % 2 == 0 { h / 2 } else { h.div_ceil(2) };
        let mut gray = Gray::new(h, weight);
        gray.with_column(al-1);

        let mut result = vec![vec![0u8; h]; h];
        
        let mut tmp = vec![0u8; h];
        for &col in &s(l-1) {
            tmp[col] ^= 1;
        }
    
        for i in (0..al - 1).rev() {
            let j = gray.previous_delta();
            // Add the row to result[j[0]] and result[j[1]]
            BinaryVector::add_inplace(&mut result[j[0]], &tmp);
            BinaryVector::add_inplace(&mut result[j[1]], &tmp);
            if i >= a {
                for &col in &s(i-a) {
                    tmp[col] ^= 1;
                }
            }
        }
        let j = gray.current_column_positions();
        for i in j {
            BinaryVector::add_inplace(&mut result[i], &tmp);
        }
        result
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use fountain_engine::traits::HDPC;
    use fountain_engine::types::CodeParams;
    use fountain_engine::DataManager;
    use fountain_utility::VecDataOperater;
    use fountain_engine::types::SolverType;
    use crate::validation::*;

    #[test]
    fn test_hdpc_10_mul_data() {
        let h = 10;
        let k = 10;
        let l = 5;
        let kl = k+l;
        let vec_length = kl;
        let params = CodeParams::new(k, k, l, h);

        //let matrix_x = gen_random_matrix(kl, vec_length, 1);
        let mut matrix_x = vec![vec![0u8; vec_length]; kl];
        for i in 0..kl {
            matrix_x[i][i] = 1;
        }

        let mut data_manager = DataManager::new_with_operator(Box::new(VecDataOperater::new(vec_length)));
        data_manager.config_from(params.clone(), SolverType::OrdEnc);

        let x_ids = (0..kl).collect::<Vec<usize>>();
        for i in x_ids.clone() {
            data_manager.insert_data_vector(i, &matrix_x[i]);
        }
        let mut y_ids = Vec::with_capacity(h);
        for _ in 0..h {
            y_ids.push(data_manager.temp_data_id());
        }
        println!("x_ids: {x_ids:?}");
        println!("y_ids: {y_ids:?}");

        let hdpc = R10HDPC::new();
        hdpc.mul_data(&mut data_manager, &params, &x_ids, &y_ids);

        for i in 0..h {
            println!("y_ids[{i}]: {:?}", data_manager.get_data_vector(y_ids[i]));
        }
    }

    fn binomial(n: usize, k: usize) -> u64 {
        if k > n {
            return 0;
        }
        if k == 0 || k == n {
            return 1;
        }
        let k = k.min(n - k); // Take advantage of symmetry
        let mut result: u64 = 1;
        for i in 0..k {
            result = result * (n - i) as u64 / (i + 1) as u64;
        }
        result
    }
    // return the smallest h such that C(h, h/2) >= k
    fn get_h(k: usize) -> usize {
        let mut h = 1;
        while binomial(h, h/2) < k as u64 {
            h += 1;
        }
        h
    }

    #[test]
    fn test_cross_validate_hdpc() {
        for k in 21..250 {
            let l = k / 10;
            let h = get_h(k+l);
            let params = CodeParams::new(k, k, l, h);
            let hdpc = R10HDPC::new();
            //validate_hdpc_mul_sparse(&hdpc, &params);
            cross_validate_hdpc(&hdpc, &params);
        }
    }
}
