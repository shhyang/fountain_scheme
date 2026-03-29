// Copyright (c) 2026 Shenghao Yang.
// All rights reserved.

use fountain_engine::traits::LDPC;

/// Validates that an LDPC implementation is consistent.
///
/// This function verifies that the LDPC implementation is consistent by checking:
/// - For each row in the active matrix, the columns returned by `active_row` 
///   contain the row index when queried via `active_column`
/// - For each row in the inactive matrix, the columns returned by `inactive_row`
///   contain the row index when queried via `inactive_column`
/// # Panics 
/// Panics if the LDPC implementation is inconsistent, with a detailed error message indicating which row/column pair failed validation.
pub fn cross_validate_ldpc<T: LDPC>(ldpc: &T, l: usize) {
    // Validate active matrix consistency
    for i in 0..l {
        let adj_checks = ldpc.active_row(i);
        for j in adj_checks {
            let adj_vars = ldpc.active_column(j);
            assert!(
                adj_vars.contains(&i),
                "LDPC validation failed: `active_row`({i}) returned column {j}, but `active_column`({j}) does not contain row {i}"
            );
        }
    }
    
    // Validate inactive matrix consistency
    for i in 0..l {
        let adj_checks = ldpc.inactive_row(i);
        for j in adj_checks {
            let adj_vars = ldpc.inactive_column(j);
            assert!(
                adj_vars.contains(&i),
                "LDPC validation failed: `inactive_row`({i}) returned column {j}, but `inactive_column`({j}) does not contain row {i}"
            );
        }
    }
}
