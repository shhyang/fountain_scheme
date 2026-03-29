// Copyright (c) 2025 Shenghao Yang. All rights reserved.
// Licensed under the MIT License. See LICENSE-MIT for details.

use rand::Rng;
use rand::prelude::SliceRandom;
use rand::rngs::StdRng;

//pub use super::degree_distribution::DegreeDistribution;
//pub use super::random::RandomDegreeSet;
//pub use super::triangular::TriangularDegreeSet;
//pub use super::random_triangular::RandomTriangularDegreeSet;
//pub use super::super::parameters::CodeParams;
//pub use super::super::traits::DegreeSet;

/// Converts a probability density function into a cumulative distribution function.
pub fn pdf_to_cdf(pdf: &Vec<f64>) -> Vec<f64> {
    let mut cdf = Vec::with_capacity(pdf.len());
    let mut acc = 0.0;
    for p in pdf.iter() {
        acc += p;
        cdf.push(acc);
    }
    // Ensure the last value is exactly 1.0 to handle floating point precision
    if !cdf.is_empty() {
        let last_idx = cdf.len() - 1;
        cdf[last_idx] = 1.0;
    }
    cdf
}

/// Truncates `pdf` so its length is at most `k + 1`, folding the tail mass into the last entry.
pub fn shorten_pdf(pdf: &mut Vec<f64>, k: usize) {
    if k == 0 {
        pdf.clear();
        return;
    }
    if pdf.len() > k + 1 {
        // sum the tail of the pdf from the (k+1)-th element to the end
        let tail = pdf.iter().skip(k).sum::<f64>();
        pdf.truncate(k + 1);
        // make the last element of the pdf the tail
        pdf[k] = tail;
    }
}

/// Uniformly samples `degree` distinct indices from `0..a` using Fisher-Yates shuffle.
pub fn sampling_indices_with_degree(
    degree: usize,
    a: usize,
    seeded_rng: &mut StdRng,
) -> Vec<usize> {
    let mut indices: Vec<usize> = (0..a).collect();
    indices.as_mut_slice().shuffle(seeded_rng);
    indices.truncate(degree);
    indices
}

/// Inverse-CDF sampler: draws a random degree from a precomputed CDF.
pub struct DegreeSampling {
    cdf: Vec<f64>,
}

impl DegreeSampling {
    /// Creates a sampler from a CDF (or PDF that will be converted).
    pub fn new(pdf: Vec<f64>) -> Self {
        let cdf = pdf_to_cdf(&pdf);
        Self { cdf }
    }

    /// Draws a single degree sample using inverse-CDF lookup.
    pub fn random(&self, seeded_rng: &mut StdRng) -> usize {
        let r = seeded_rng.r#gen::<f64>();
        self.cdf.iter().position(|&x| r < x).unwrap_or(1)
    }
}

#[cfg(test)]
mod tests {
    //use super::*;
    //use super::super::degree_distribution::DegreeDistribution;
    //use super::super::random::RandomDegreeSet;
    //use super::super::triangular::TriangularDegreeSet;
    //use crate::fountain_code::parameters::{CodeParams, SolverType};
    //use crate::fountain_code::precode::PrecodeType;
    use crate::degree_sets::SolitonDistribution;
    use crate::degree_sets::DegreeSetGenerator;
    use crate::degree_sets::random::RandomDegreeSet;
    use crate::degree_sets::triangular::TriangularDegreeSet;
    use crate::parameters::ParamsConfig;
    #[test]
    fn test_ideal_soliton_pdf_sum_to_one() {
        let k = 10;
        let dist = SolitonDistribution::IdealSoliton { dmax: k };
        let pdf = dist.pdf();
        let sum: f64 = pdf.iter().sum();
        assert!(
            (sum - 1.0).abs() < 1e-8,
            "PDF does not sum to 1, got {}",
            sum
        );
    }

    #[test]
    fn test_ideal_soliton_pdf_is_correct() {
        let k = 10;
        let dist = SolitonDistribution::IdealSoliton { dmax: k };
        let pdf = dist.pdf();
        let expected_pdf = vec![
            0.0,
            0.1,
            0.5,
            1.0 / 6.0,
            1.0 / 12.0,
            1.0 / 20.0,
            1.0 / 30.0,
            1.0 / 42.0,
            1.0 / 56.0,
            1.0 / 72.0,
            1.0 / 90.0,
        ];
        assert_eq!(pdf, expected_pdf);
    }

    #[test]
    fn test_robust_soliton_pdf_sum_to_one() {
        let k = 10;
        let dist = SolitonDistribution::RobustSoliton {
            dmax: k,
            c: 0.1,
            delta: 0.05,
        };
        let pdf = dist.pdf();
        println!("PDF: {:?}", pdf);
        let sum: f64 = pdf.iter().sum();
        assert!(
            (sum - 1.0).abs() < 1e-8,
            "PDF does not sum to 1, got {}",
            sum
        );
    }

    #[test]
    fn test_ideal_random_degree_set_degree_range() {
        let k = 100;
        let dmax = 120;
        let dist = SolitonDistribution::IdealSoliton { dmax };
        let params = ParamsConfig::LT.generate(k);
        let mut assoc = RandomDegreeSet::new(&dist, &params);
        for coded_id in 0..100 {
            let (active_indices, _inactive_indices) = assoc.degree_set(coded_id);
            assert!(!active_indices.is_empty() && active_indices.len() <= k);
        }
    }

    #[test]
    fn test_robust_random_degree_set_degree_range() {
        let k = 100;
        let dist = SolitonDistribution::RobustSoliton {
            dmax: k,
            c: 0.1,
            delta: 0.5,
        };
        let params = ParamsConfig::LT.generate(k);
        let mut assoc = RandomDegreeSet::new(&dist, &params);
        for coded_id in 0..100 {
            let (active_indices, _inactive_indices) = assoc.degree_set(coded_id);
            assert!(!active_indices.is_empty() && active_indices.len() <= k);
        }
    }

    #[test]
    fn test_ideal_uniform_sampling_of_indices() {
        let k = 20;
        let ideal = SolitonDistribution::IdealSoliton { dmax: k };
        let params = ParamsConfig::LT.generate(k);
        let mut assoc = RandomDegreeSet::new(&ideal, &params);

        let num_trials = 10_000;
        let mut counts = vec![0usize; k];

        for coded_id in 0..num_trials {
            let set = assoc.degree_set(coded_id);
            for &idx in &set.0 {
                counts[idx] += 1;
            }
        }

        let avg = counts.iter().sum::<usize>() as f64 / k as f64;
        let max_dev = counts
            .iter()
            .map(|&c| ((c as f64 - avg).abs() / avg))
            .fold(0.0, f64::max);

        println!("Counts: {:?}", counts);
        println!("Average: {}", avg);
        println!("Max relative deviation: {:.2}%", max_dev * 100.0);

        // Allow up to 10% deviation from the mean (statistical fluctuations)
        assert!(max_dev < 0.10, "Sampling is not uniform enough");
    }

    #[test]
    fn test_robust_uniform_sampling_of_indices() {
        let k = 20;
        let dist = SolitonDistribution::RobustSoliton {
            dmax: k,
            c: 0.1,
            delta: 0.5,
        };
        let params = ParamsConfig::LT.generate(k);
        let mut assoc = RandomDegreeSet::new(&dist, &params);
        let num_trials = 10_000;
        let mut counts = vec![0usize; k];

        for coded_id in 0..num_trials {
            let set = assoc.degree_set(coded_id);
            for &idx in &set.0 {
                counts[idx] += 1;
            }
        }

        let avg = counts.iter().sum::<usize>() as f64 / k as f64;
        let max_dev = counts
            .iter()
            .map(|&c| ((c as f64 - avg).abs() / avg))
            .fold(0.0, f64::max);

        println!("Counts: {:?}", counts);
        println!("Average: {}", avg);
        println!("Max relative deviation: {:.2}%", max_dev * 100.0);
    }

    #[test]
    fn test_triangular_degree_set_degree_range() {
        let k = 100;
        let dmax = 120;
        let dist = SolitonDistribution::IdealSoliton { dmax };
        let params = ParamsConfig::LT.generate(k);
        let mut assoc = TriangularDegreeSet::new(dist, params);
        for coded_id in 0..100 {
            let (active_indices, _inactive_indices) = assoc.degree_set(coded_id);
            assert!(
                !active_indices.is_empty() && active_indices.iter().max().unwrap() == &coded_id
            );
        }
    }
}
