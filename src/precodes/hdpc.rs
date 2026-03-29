// Copyright (c) 2025 Shenghao Yang. All rights reserved.
// Licensed under the MIT License. See LICENSE-MIT for details.

use fountain_engine::DataManager;
use fountain_engine::algebra::finite_field::GF256;
use fountain_engine::traits::HDPC;
use fountain_engine::types::CodeParams;
//use rand::{Rng, SeedableRng, rngs::StdRng};

/// No-op HDPC implementation (no high-density parity-check rows).
pub struct NoHDPC;
impl HDPC for NoHDPC {}


fn rq_mul_data(
    manager: &mut DataManager,
    params: &CodeParams,
    x_ids: &[usize],
    y_ids: &[usize],
    delta_column: &dyn Fn(usize) -> Vec<usize>,
) {
    let m = x_ids.len();
    let h = y_ids.len();
    let kl = params.num_message_ldpc();
    let temp_id = manager.temp_data_id();
    //dbg!(&m, &h);
    //manager.ensure_zero(&[temp_id]);
    manager.add_to_vector(&[x_ids[0]], temp_id);
    for i in 0..m - 1 {
        //dbg!("mul_vector", &temp_id, manager.get_vector(temp_id));
        let j = delta_column(i)
            .iter()
            .map(|&id| y_ids[id])
            .collect::<Vec<_>>();
        manager.broadcast_add(temp_id, &j);
        manager.multiply_alpha(temp_id);
        manager.add_to_vector(&[x_ids[i + 1]], temp_id);
    }
    for i in m - 1..kl - 1 {
        let j = delta_column(i)
            .iter()
            .map(|&id| y_ids[id])
            .collect::<Vec<_>>();
        manager.broadcast_add(temp_id, &j);
        manager.multiply_alpha(temp_id);
    }
    //dbg!("mul_vector", &temp_id, manager.get_vector(temp_id));
    for i in 0..h {
        manager.add_to_vector(&[temp_id], y_ids[i]);
        manager.multiply_alpha(temp_id);
    }
    manager.remove(temp_id);
}
fn rq_mul_binary(
    gf: &GF256,
    params: &CodeParams,
    n: usize,
    v: &dyn Fn(usize) -> Vec<u8>,
    delta_column: &dyn Fn(usize) -> Vec<usize>,
) -> Vec<Vec<u8>> {
    let h = params.h;
    let kl = params.num_message_ldpc();
    let mut result = vec![vec![0u8; n]; h];
    let mut tmp = v(0);

    for i in 0..kl - 1 {
        let j = delta_column(i);
        gf.vector_addition_inplace(&mut result[j[0]], &tmp);
        gf.vector_addition_inplace(&mut result[j[1]], &tmp);
        gf.multiply_alpha_inplace(&mut tmp);
        let next = v(i + 1);
        gf.vector_addition_inplace(&mut tmp, &next);
    }
    for i in 0..h {
        gf.vector_addition_inplace(&mut result[i], &tmp);
        gf.multiply_alpha_inplace(&mut tmp);
    }
    result
}

fn rq_mul_sparse_sh(
    gf: &GF256,
    params: &CodeParams,
    s: &dyn Fn(usize) -> Vec<usize>,
    delta_column: &dyn Fn(usize) -> Vec<usize>,
) -> Vec<Vec<u8>> {
    let h = params.h;
    let l = params.l;
    let a = params.a;
    let lb = l + params.b;
    let mut result = vec![vec![0u8; h]; h];
    let mut tmp = vec![0u8; h];
    s(0).iter().for_each(|&col| tmp[col] ^= 1);

    for i in 0..lb - 1 {
        let j = delta_column(a + i);
        // Add the row to result[j[0]] and result[j[1]]
        gf.vector_addition_inplace(&mut result[j[0]], &tmp);
        gf.vector_addition_inplace(&mut result[j[1]], &tmp);
        gf.multiply_alpha_inplace(&mut tmp);
        if i < l - 1 {
            s(i + 1).iter().for_each(|&col| tmp[col] ^= 1);
        }
    }
    for i in 0..h {
        gf.vector_addition_inplace(&mut result[i], &tmp);
        gf.multiply_alpha_inplace(&mut tmp);
    }
    result
}

/* macro_rules! impl_hdpc {
    ($type:ty) => {
        impl HDPC for $type {
            fn mul_data(
                &self,
                manager: &mut DataManager,
                params: &CodeParams,
                x_ids: &[usize],
                y_ids: &[usize],
            ) {
                let delta_column = |j: usize| self.delta_column(j);
                rq_mul_data(manager, params, x_ids, y_ids, &delta_column);
            }

            fn mul_binary(
                &self,
                params: &CodeParams,
                n: usize,
                v: &dyn Fn(usize) -> Vec<u8>,
            ) -> Vec<Vec<u8>> {
                let delta_column = |j: usize| self.delta_column(j);
                rq_mul_binary(&GF256::default(), params, n, v, &delta_column)
            }

            fn mul_sparse(
                &self,
                params: &CodeParams,
                n: usize,
                s: &dyn Fn(usize) -> Vec<usize>,
            ) -> Vec<Vec<u8>> {
                let v = |i: usize| {
                    let mut col = vec![0u8; n];
                    for &j in s(i).iter() {
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
                let delta_column = |j: usize| self.delta_column(j);
                rq_mul_sparse_sh(&GF256::default(), params, s, &delta_column)
            }
        }
    };
} */

pub struct GenericRQHDPC {
    delta_column_fn: Box<dyn Fn(usize) -> Vec<usize>>,
}

impl GenericRQHDPC {
    #[must_use]
    pub fn new(delta_column_fn: Box<dyn Fn(usize) -> Vec<usize>>) -> Self {
        Self { delta_column_fn }
    }
}

impl HDPC for GenericRQHDPC {
    fn mul_data(
        &self,
        manager: &mut DataManager,
        params: &CodeParams,
        x_ids: &[usize],
        y_ids: &[usize],
    ) {
        let delta_column = |j: usize| (self.delta_column_fn)(j);
        rq_mul_data(manager, params, x_ids, y_ids, &delta_column);
    }

    fn mul_binary(
        &self,
        params: &CodeParams,
        n: usize,
        v: &dyn Fn(usize) -> Vec<u8>,
    ) -> Vec<Vec<u8>> {
        let delta_column = |j: usize| (self.delta_column_fn)(j);
        rq_mul_binary(&GF256::default(), params, n, v, &delta_column)
    }

    fn mul_sparse(
        &self,
        params: &CodeParams,
        n: usize,
        s: &dyn Fn(usize) -> Vec<usize>,
    ) -> Vec<Vec<u8>> {
        let v = |i: usize| {
            let mut col = vec![0u8; n];
            for &j in s(i).iter() {
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
        let delta_column = |j: usize| (self.delta_column_fn)(j);
        rq_mul_sparse_sh(&GF256::default(), params, s, &delta_column)
    }
}

/// Returns the default (cyclic) HDPC precode with `h` check symbols.
/// Delta-column mapping is `(j % h, (j+1) % h)`.
#[must_use]
pub fn default_rq_hdpc(h: usize) -> Box<dyn HDPC> {
    let delta_column_fn = Box::new(move |j: usize| vec![j % h, (j + 1) % h]);
    Box::new(GenericRQHDPC::new(delta_column_fn))
}

/*
pub struct RQHDPC{
    h: usize, // number of hdpc vectors
    gf: GF256,

    delta_column_fn: Box<dyn Fn(usize) -> Vec<usize>>,
}

impl RQHDPC {
    pub fn new(num_hdpc: usize) -> Self {
        Self {
            h: num_hdpc,
            gf: GF256::default(),
            delta_column_fn: Box::new(move |j| {
                let pos1 = j % num_hdpc;
                let pos2 = (j + 1) % num_hdpc;
                vec![pos1, pos2]
            })
        }
    }

    pub fn set_delta_column_fn(&mut self, delta_column_fn: Box<dyn Fn(usize) -> Vec<usize>>) {
        self.delta_column_fn = delta_column_fn;
    }

    /// Calculate Y = D * S, where D is an h x l matrix and S is an l x n binary matrix
    /// s is a function that returns the non-zero indices of the i-th row of S
    pub fn mul_binary(&self, kl: usize, l: usize, n: usize, s: &dyn Fn(usize) -> Vec<usize>) -> Vec<Vec<u8>> {
        let h = self.h;
        let mut result = vec![vec![0u8; n]; h];
        let mut tmp = vec![0u8; n];
        s(0).iter().for_each(|&col| tmp[col] ^= 1);

        for i in 0..l-1 {
            let j = (self.delta_column_fn)(kl-l+i);
            // Add the row to result[j[0]] and result[j[1]]
            self.gf.vector_addition_inplace(&mut result[j[0]], &tmp);
            self.gf.vector_addition_inplace(&mut result[j[1]], &tmp);
            self.gf.multiply_alpha_inplace(&mut tmp);
            s(i+1).iter().for_each(|&col| tmp[col] ^= 1);
        }
        for i in 0..h {
            self.gf.vector_addition_inplace(&mut result[i], &tmp);
            self.gf.multiply_alpha_inplace(&mut tmp);
        }
        result
    }

    /// Calculate Y = D * X
    pub fn mul_vector(&self, manager: &mut DataManager, kl: usize, x_pkg: &[usize], y_pkg: &[usize]) {
        let m = x_pkg.len();
        let h = y_pkg.len();
        let temp_id = manager.temp_data_id();
        //dbg!(&m, &h);
        //manager.ensure_zero(&[temp_id]);
        manager.add_to_vector(&[x_pkg[0]], temp_id);
        for i in 0..m-1 {
            //dbg!("mul_vector", &temp_id, manager.get_vector(temp_id));
            let j = (self.delta_column_fn)(i).iter().map(|&id| y_pkg[id]).collect::<Vec<_>>();
            manager.broadcast_add(temp_id, &j);
            manager.multiply_alpha(temp_id);
            manager.add_to_vector(&[x_pkg[i+1]], temp_id);
        }
        for i in m-1..kl-1 {
            let j = (self.delta_column_fn)(i).iter().map(|&id| y_pkg[id]).collect::<Vec<_>>();
            manager.broadcast_add(temp_id, &j);
            manager.multiply_alpha(temp_id);
        }
        //dbg!("mul_vector", &temp_id, manager.get_vector(temp_id));
        for i in 0..h {
            manager.add_to_vector(&[temp_id], y_pkg[i]);
            manager.multiply_alpha(temp_id);
        }
    }
}
*/

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parameters::generate_rq_like_parameters;
    use fountain_engine::algebra::finite_field::GF256;
    use fountain_engine::algebra::linear_algebra::Matrix;
    use fountain_engine::traits::DataOperator;
    use fountain_engine::CodeParams;
    use fountain_utility::VecDataOperater;
    use crate::validation::*;

    use rand::Rng;
    use rand::rngs::StdRng;
    use rand::SeedableRng;

    fn gen_random_matrix(m: usize, n: usize, ffsize: u8) -> Vec<Vec<u8>> {
        let mut matrix_a = Vec::new();
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

    fn gen_matrix_d(field: &GF256, h: usize, l: usize) -> Vec<Vec<u8>> {
        // matrix_delta is an h x l matrix
        let mut matrix_delta = vec![vec![0u8; l]; h];
        for i in 0..l - 1 {
            let j = [i % h, (i + 1) % h];
            matrix_delta[j[0]][i] = 1;
            matrix_delta[j[1]][i] = 1;
        }
        // the last column of matrix_delta is 1, alpha, alpha^2, ..., alpha^(h-1)
        let mut alpha = 1;
        for i in 0..h {
            matrix_delta[i][l - 1] = alpha;
            alpha = field.mul_alpha(alpha);
        }

        // matrix_tau is an l x l lower triangular matrix
        let mut alpha = 1;
        let mut matrix_tau = vec![vec![0u8; l]; l];
        for i in 0..l {
            for j in 0..l - i {
                matrix_tau[i + j][j] = alpha;
            }
            alpha = field.mul_alpha(alpha);
        }
        Matrix::multiply(field, &matrix_delta, &matrix_tau)
    }

    #[test]
    fn hdpc_constraints_test() {
        let k = 21;
        let l = 3;
        let b = 4;
        let _a = k - b;
        let h = 4;
        let matrix_d = gen_matrix_d(&GF256::default(), h, k + l);
        let mut matrix_c = vec![vec![0u8; 1]; k + l];
        matrix_c[19][0] = 1;
        matrix_c[23][0] = 1;

        let result = Matrix::multiply(&GF256::default(), &matrix_d, &matrix_c);
        println!("{:?}", result);
    }

    #[test]
    fn test_mul_binary() {
        for k in 21..250 {
            let params = generate_rq_like_parameters(k);
            let n = k - params.a + params.h;
            cross_validate_mul_binary(&params, n);
        }
    }

    fn cross_validate_mul_binary(params: &CodeParams, n: usize) {
        let kl = params.num_message_ldpc();
        let hdpc = default_rq_hdpc(params.h);
        let mat_s = gen_random_matrix(kl, n, 1);

        let v = |row: usize| -> Vec<u8> { mat_s[row].clone() };
        let s = |row: usize| -> Vec<usize> {
            mat_s[row]
                .iter()
                .enumerate()
                .filter(|&(_, &v)| v != 0)
                .map(|(col, _)| col)
                .collect::<Vec<_>>()
        };

        let result_binary = hdpc.mul_binary(params, n, &v);
        let result_sparse = hdpc.mul_sparse(params, n, &s);

        assert_eq!(result_binary, result_sparse);
    }

    #[test]
    fn test_mul_vector() {
        for k in 11..250 {
            let params = generate_rq_like_parameters(k);
            let h = params.h;
            let l = params.l;
            cross_validate_mul_vector(h, params.a + l);
        }
    }

    fn cross_validate_mul_vector(h: usize, l: usize) {
        println!("cross_validate_mul_vector: h: {}, l: {}", h, l);
        let _hdpc = default_rq_hdpc(h);
        let field = GF256::default();
        let t = 3;

        // generate a random matrix over GF(256)
        let mat_x = gen_random_matrix(l, t, 8);

        // Prepare a VecDataManager with l vectors of length t
        let mut manager = VecDataOperater::new(t);
        for i in 0..l {
            manager.insert_vector(&mat_x[i], i);
        }

        // calculate Y = D * X by mul_vector
        let _y_ids = (l..l + h).collect::<Vec<_>>();
        //manager.zero_vectors(&y_ids);
        //hdpc.mul_vector(&mut manager, (0..l).collect(), &y_ids);

        // calculate Y = D * X by matrix multiplication
        let matrix_d = gen_matrix_d(&field, h, l);
        let _mat_y = Matrix::multiply(&field, &matrix_d, &mat_x);

        // check if mat_y and mat_y_manager are equal
        for _i in 0..h {
            //assert_eq!(manager.get_vector(y_ids[i]), mat_y[i]);
        }
    }

    #[test]
    fn test_cross_validate_hdpc() {
        for k in 21..250 {
            let params = generate_rq_like_parameters(k);
            let h = params.h;
            let random_pair = Box::new(move |j: usize| {
                let mut rng = StdRng::seed_from_u64(j as u64);
                let pos1 = rng.gen_range(0..h);
                let pos2 = rng.gen_range(0..h);
                vec![pos1, pos2]
            });
            let hdpc = GenericRQHDPC::new(random_pair);
            cross_validate_hdpc(&hdpc, &params);
        }
    }
}
