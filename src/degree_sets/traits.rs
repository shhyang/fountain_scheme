// Copyright (c) 2025 Shenghao Yang. All rights reserved.
// Licensed under the MIT License. See LICENSE-MIT for details.

/// Trait for generating degree sets for fountain codes
///
/// This trait defines the interface for degree set generators used in
/// fountain coding algorithms. Degree sets determine how many source
/// symbols each coded symbol should be connected to.
pub trait DegreeSetGenerator {
    /// Generate a degree set for a given coded symbol ID
    ///
    /// # Arguments
    ///
    /// * `coded_id` - The ID of the coded symbol
    ///
    /// # Returns
    ///
    /// A tuple containing:
    /// - `Vec<usize>`: Active indices (source symbols to connect to)
    /// - `Vec<usize>`: Inactive indices (for future use)
    fn degree_set(&mut self, coded_id: usize) -> (Vec<usize>, Vec<usize>);

    /// Get the name of this degree set generator
    ///
    /// # Returns
    ///
    /// A string slice containing the name of the generator
    fn name(&self) -> &str;
}
