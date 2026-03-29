// Copyright (c) 2025 Shenghao Yang. All rights reserved.
// Licensed under the MIT License. See LICENSE-MIT for details.

use super::traits::DegreeSetGenerator;
use super::{SolitonDistribution, sampling::shorten_pdf};
use fountain_engine::CodeParams;

/// Degree set using triangular embedding.
pub struct TriangularDegreeSet {
    pdf: Vec<f64>,
    params: CodeParams, // range of indices for sampling
    degree_sampled: Vec<f64>,
    indices_sampled: Vec<usize>,
}

impl TriangularDegreeSet {
    /// Creates a new triangular degree set from the distribution and code parameters.
    pub fn new(degree_distribution: SolitonDistribution, params: CodeParams) -> Self {
        if params.k == 0 {
            panic!("Source count must be greater than 0");
        }

        let mut pdf = degree_distribution.pdf();
        shorten_pdf(&mut pdf, params.a);
        let dmax = pdf.len();

        Self {
            pdf,
            params: params.clone(),
            degree_sampled: vec![0.0; dmax],
            indices_sampled: vec![0; params.a],
        }
    }

    /// Generates a lower-triangular index set: the returned set always includes `coded_id`
    /// as its maximum element, ensuring the first `a` coded vectors form a lower-triangular matrix.
    pub fn triangular(&mut self, coded_id: usize) -> Vec<usize> {
        // add pdf[i] to degree_sampled[i] for all i >= 0
        self.degree_sampled
            .iter_mut()
            .enumerate()
            .for_each(|(i, p)| *p += self.pdf[i]);

        if coded_id == 0 {
            self.degree_sampled[1] -= 1.0;
            self.indices_sampled[0] += 1;
            return vec![0];
        }
        if coded_id == 1 {
            self.degree_sampled[2] -= 1.0;
            self.indices_sampled[0] += 1;
            self.indices_sampled[1] += 1;
            return vec![0, 1];
        }
        // degree is the index d such that degree_sampled[d] is the maximum.
        let degree = self
            .degree_sampled
            .iter()
            .enumerate()
            .max_by(|(_, p1), (_, p2)| p1.partial_cmp(p2).unwrap())
            .map(|(i, _)| i)
            .unwrap();
        self.degree_sampled[degree] -= 1.0;
        let mut indices = self.indices_sampled.clone();
        if coded_id >= 2 && coded_id < self.params.a {
            indices.truncate(coded_id);
            indices.as_mut_slice().sort();
            indices.truncate(degree - 1);
            indices.push(coded_id);
        } else {
            indices.as_mut_slice().sort();
            indices.truncate(degree);
        }
        indices
            .iter_mut()
            .for_each(|i| self.indices_sampled[*i] += 1);
        indices
    }
}

impl DegreeSetGenerator for TriangularDegreeSet {
    fn degree_set(&mut self, coded_id: usize) -> (Vec<usize>, Vec<usize>) {
        let active_indices = self.triangular(coded_id);
        let inactive_indices = vec![];
        (active_indices, inactive_indices)
    }
    fn name(&self) -> &str {
        "Triangular"
    }
}
