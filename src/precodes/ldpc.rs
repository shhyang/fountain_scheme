// Copyright (c) 2025 Shenghao Yang. All rights reserved.
// Licensed under the MIT License. See LICENSE-MIT for details.

use fountain_engine::CodeParams;
use fountain_engine::traits::LDPC;

/// No-op LDPC implementation (identity — no parity-check rows added).
pub struct NoLDPC;
impl LDPC for NoLDPC {}

/// LDPC precode with reversed column order, suitable for systematic encoding.
#[derive(Debug, Clone)]
pub struct ReversedLDPC {
    l: usize,
    a: usize,
    /// Total inactive count: `b + h`.
    n: usize,
}

impl ReversedLDPC {
    /// Creates a new reversed-order LDPC precode from the given code parameters.
    #[must_use]
    pub fn new(params: &CodeParams) -> Self {
        Self {
            l: params.l,
            a: params.a,
            n: params.b + params.h,
        }
    }
}

impl LDPC for ReversedLDPC {
    /// Non-zero entries for a column of `S_a`, which is an l x a matrix
    #[inline]
    fn active_column(&self, var_col: usize) -> Vec<usize> {
        let l = self.l;
        let rev_col = self.a - var_col - 1;
        let which_block = rev_col / l;
        let pos_in_block = rev_col % l;
        let mut base_positions = vec![0, (which_block + 1) % l, (2 * which_block + 2) % l];
        base_positions.sort_unstable();
        base_positions.dedup();
        base_positions
            .iter()
            .map(|&base_pos| (base_pos + pos_in_block) % l)
            .collect()
    }

    /// Non-zero entries for a row of `S_a`, which is an l x a matrix
    #[inline]
    fn active_row(&self, check_row: usize) -> Vec<usize> {
        let l = self.l;
        let mut var_cols = Vec::new();
        let num_blocks = self.a / l;
        let extra_cols = self.a - num_blocks * l;
        for i in 0..num_blocks {
            let check_ids = self.active_column(self.a + check_row - (i + 1) * l);
            var_cols.extend(
                check_ids
                    .iter()
                    .map(|&id| id + extra_cols + (num_blocks - i - 1) * l),
            );
        }
        if extra_cols > 0 {
            //let check_ids = self.active_column(extra_cols + check_row - l); // need to use isize
            //var_cols.extend(check_ids.iter().map(|&id| id as isize + extra_cols as isize - l as isize).filter(|&id| id >= 0).map(|id| id as usize));
            let rev_col = self.a - 1 + l - extra_cols - check_row;
            let which_block = rev_col / l;
            let pos_in_block = rev_col % l;
            let mut base_positions = vec![0, (which_block + 1) % l, (2 * which_block + 2) % l];
            base_positions.sort_unstable();
            base_positions.dedup();
            var_cols.extend(
                base_positions
                    .iter()
                    .map(|&base_pos| (base_pos + pos_in_block) % l)
                    .filter(|&id| id >= l - extra_cols)
                    .map(|id| id + extra_cols - l),
            );
        }
        var_cols
    }

    /// Non-zero entries for a row of `S_i`, which is an l x (b+h) matrix
    #[inline]
    fn inactive_row(&self, check_row: usize) -> Vec<usize> {
        let n = self.n;
        if n == 0 {
            return vec![];
        }
        vec![n - 1 - ((check_row + 1) % n), n - 1 - (check_row % n)]
    }
    /// Non-zero entries for a column of `S_i`, which is an l x (b+h) matrix
    #[inline]
    fn inactive_column(&self, var_col: usize) -> Vec<usize> {
        let l = self.l;
        let n = self.n;
        if n == 0 {
            return vec![];
        }
        let num_blocks = l / n;
        let check_ids = self.inactive_row(var_col);
        let mut check_rows = Vec::new();
        for i in 0..num_blocks {
            check_rows.extend(check_ids.iter().map(|&id| id + i * n));
        }
        check_rows.extend(
            check_ids
                .iter()
                .map(|&id| id + num_blocks * n)
                .filter(|&id| id < l),
        );
        check_rows
    }
}

/// Standard RaptorQ LDPC precode (forward column order), delegates to [`ReversedLDPC`] with index reversal.
#[derive(Debug)]
pub struct RQLDPC {
    reversed_ldpc: ReversedLDPC,
}

impl RQLDPC {
    /// Creates a new standard RaptorQ LDPC precode from the given code parameters.
    #[must_use]
    pub fn new(params: &CodeParams) -> Self {
        Self {
            reversed_ldpc: ReversedLDPC::new(params),
        }
    }
}

impl LDPC for RQLDPC {
    /// Non-zero entries for a column of `S_a`
    #[inline]
    fn active_column(&self, var_col: usize) -> Vec<usize> {
        if self.reversed_ldpc.a == 0 {
            return vec![];
        }
        self.reversed_ldpc
            .active_column(self.reversed_ldpc.a - 1 - var_col)
    }

    /// Non-zero entries for a row of `S_i`
    #[inline]
    fn active_row(&self, check_row: usize) -> Vec<usize> {
        if self.reversed_ldpc.a == 0 {
            return vec![];
        }
        let aminus1 = self.reversed_ldpc.a - 1;
        self.reversed_ldpc
            .active_row(check_row)
            .iter()
            .map(|&id| aminus1 - id)
            .collect()
    }

    /// Non-zero entries for a row of `S_i`
    #[inline]
    fn inactive_row(&self, check_row: usize) -> Vec<usize> {
        let n = self.reversed_ldpc.n;
        if n == 0 {
            return vec![];
        }
        vec![check_row % n, (check_row + 1) % n]
    }
    /// Non-zero entries for a column of `S_b`
    #[inline]
    fn inactive_column(&self, var_col: usize) -> Vec<usize> {
        if self.reversed_ldpc.n == 0 {
            return vec![];
        }
        self.reversed_ldpc
            .inactive_column(self.reversed_ldpc.n - 1 - var_col)
    }
}

/// Raptor-10 LDPC variant (same structure as [`RQLDPC`] but without inactive-column methods).
pub struct R10LDPC {
    reversed_ldpc: ReversedLDPC,
}

impl R10LDPC {
    /// Creates a new Raptor-10 LDPC precode from the given code parameters.
    #[must_use]
    pub fn new(params: &CodeParams) -> Self {
        Self {
            reversed_ldpc: ReversedLDPC::new(params),
        }
    }
}

impl LDPC for R10LDPC {
    /// Non-zero entries for a column of `S_a`
    #[inline]
    fn active_column(&self, var_col: usize) -> Vec<usize> {
        self.reversed_ldpc
            .active_column(self.reversed_ldpc.a - 1 - var_col)
    }

    /// Non-zero entries for a row of `S_i`
    #[inline]
    fn active_row(&self, check_row: usize) -> Vec<usize> {
        let aminus1 = self.reversed_ldpc.a - 1;
        self.reversed_ldpc
            .active_row(check_row)
            .iter()
            .map(|&id| aminus1 - id)
            .collect()
    }
}

/// Test the precode
#[cfg(test)]
#[allow(clippy::many_single_char_names)]
mod tests {
    use super::*;
    use crate::parameters::generate_rq_like_parameters;
    use crate::validation::cross_validate_ldpc;
    use fountain_engine::algebra::finite_field::GF256;
    use fountain_engine::algebra::linear_algebra::Matrix;
    use fountain_engine::traits::LDPC;

    fn gen_parity_check_matrix<T: LDPC>(ldpc: &T, a: usize, l: usize, n: usize) -> Vec<Vec<u8>> {
        let total_vars = a + n + l;

        let mut matrix = vec![vec![0u8; total_vars]; l];

        // Fill the matrix by checking each variable node's connections
        for var_id in 0..a {
            let adj_checks = ldpc.active_column(var_id);
            for &check_id in &adj_checks {
                matrix[check_id][var_id] = 1;
            }
        }

        for var_id in a..a + l {
            matrix[var_id - a][var_id] = 1;
        }

        for var_id in a + l..a + l + n {
            let adj_checks = ldpc.inactive_column(var_id - a - l);
            for &check_id in &adj_checks {
                matrix[check_id][var_id] = 1;
            }
        }

        matrix
    }

    #[test]
    fn ldpc_constraints_test() {
        let params = CodeParams::new(21, 17, 3, 4);
        let k = 21;
        let l = 3;
        let b = 4;
        let a = k - b;
        let h = 4;
        let ldpc = RQLDPC::new(&params);
        let matrix_h = gen_parity_check_matrix(&ldpc, a, l, b + h);
        let result = show_parity_check_matrix(&matrix_h, a, l, a + l + b + h);
        println!("{}", result);

        let mut matrix_c = vec![vec![0u8; 1]; k + l + h];
        matrix_c[19][0] = 1;
        matrix_c[23][0] = 1;
        matrix_c[24][0] = 18;
        matrix_c[25][0] = 36;
        matrix_c[26][0] = 72;
        matrix_c[27][0] = 129;
        //matrix_c[k+4][0] = 127;
        //matrix_c[k+5][0] = 42;
        //matrix_c[k+6][0] = 218;
        let result = Matrix::multiply(&GF256::default(), &matrix_h, &matrix_c);
        println!("{:?}", result);
    }

    #[test]
    fn ldpc_cross_validation_rqldpc() {
        for k in 11..150 {
            let params = generate_rq_like_parameters(k);
            let ldpc = RQLDPC::new(&params);
            cross_validate_ldpc(&ldpc, params.l);
        }
    }

    #[test]
    fn ldpc_cross_validation_reversed_ldpc() {
        for k in 11..150 {
            let params = generate_rq_like_parameters(k);
            let ldpc = ReversedLDPC::new(&params);
            cross_validate_ldpc(&ldpc, params.l);
        }
    }

    #[test]
    fn ldpc_cross_validation_r10_ldpc() {
        for k in 11..150 {
            let params = generate_rq_like_parameters(k);
            let ldpc = R10LDPC::new(&params);
            cross_validate_ldpc(&ldpc, params.l);
        }
    }
    //for i in 0..l {
    //    let adj_checks = ldpc.active_row(i);
    //    for j in adj_checks {
    //        let adj_vars = ldpc.active_column(j);
    //        assert!(adj_vars.contains(&i));
    //    }
    //}
    //for i in 0..l {
    //    let adj_checks = ldpc.inactive_row(i);
    //    dbg!(&i, &adj_checks);
    //    for j in adj_checks {
    //        let adj_vars = ldpc.inactive_column(j);
    //        dbg!(&j, &adj_vars);
    //        assert!(adj_vars.contains(&i));
    //    }
    //}

    #[test]
    fn show_parity_check_matrix_test() {
        let a = 21;
        let l = 3;
        let n = 4;
        let params = CodeParams::new(a + l + n, a, l, n);
        let ldpc = RQLDPC::new(&params);

        let matrix = gen_parity_check_matrix(&ldpc, a, l, n);

        let result = show_parity_check_matrix(&matrix, a, l, a + l + n);
        println!("{}", result);
    }

    /// Display the parity check matrix for debugging purposes
    fn show_parity_check_matrix(
        matrix: &[Vec<u8>],
        a: usize,
        l: usize,
        total_vars: usize,
    ) -> String {
        // Convert to string representation
        let mut result = String::new();
        result.push_str(&format!(
            "LDPC Parity Check Matrix ({l}×{total_vars}):\n    "
        ));

        // Print column headers
        for col in 0..total_vars {
            if col < a {
                result.push_str(&format!("A{col:<2} "));
            } else if col < a + l {
                result.push_str(&format!("L{:<2} ", col - a));
            } else {
                result.push_str(&format!("B{:<2} ", col - a - l));
            }
        }
        result.push('\n');

        // Print each row with its label and values
        for (row, row_data) in matrix.iter().enumerate() {
            result.push_str(&format!("C{row:<2} "));
            for &val in row_data {
                result.push_str(&format!(" {val}  "));
            }
            result.push('\n');
        }
        result
    }
}
