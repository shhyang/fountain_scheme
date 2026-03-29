// Copyright (c) 2025 Shenghao Yang. All rights reserved.
// Licensed under the MIT License. See LICENSE-MIT for details.

use super::degree_distribution::SolitonDistribution;
use super::examples::{DeterministicDegreeSet, FixedDegreeSet};
use super::random::RandomDegreeSet;
use super::random_triangular::RandomTriangularDegreeSet;
use super::traits::DegreeSetGenerator;
use fountain_engine::CodeParams;

/// Selects the degree-set generation strategy and its distribution parameters.
#[derive(Clone, Debug)]
pub enum DegreeSetConfig {
    /// Random sampling from the ideal soliton distribution (full range).
    RandomIdeal,
    /// Random sampling from the robust soliton distribution.
    RandomRobust {
        /// Constant that scales the spike in the robust soliton distribution.
        c: f64,
        /// Failure probability bound for the robust soliton distribution.
        delta: f64,
    },
    /// Random sampling from the ideal soliton distribution, truncated at `dmax`.
    RandomIdealTruncated {
        /// Maximum degree (capped at the active symbol count).
        dmax: usize,
    },
    /// Triangular embedding with ideal soliton distribution, truncated at `dmax`.
    TriangularIdealTruncated {
        /// Maximum degree (capped at the active symbol count).
        dmax: usize,
    },
    /// Lookup a degree-set strategy by name at runtime.
    ByName {
        /// Strategy name (e.g. `"random_ideal"`, `"deterministic"`).
        name: String,
        /// Extra numeric parameters forwarded to the strategy constructor.
        degree_params: Vec<f64>,
    },
}

impl DegreeSetConfig {
    /// Instantiates the configured degree-set generator for the given code parameters.
    pub fn create(&self, params: &CodeParams) -> Box<dyn DegreeSetGenerator> {
        match self {
            DegreeSetConfig::RandomIdeal => Box::new(RandomDegreeSet::new(
                &SolitonDistribution::IdealSoliton {
                    dmax: params.num_active(),
                },
                params,
            )),
            DegreeSetConfig::RandomRobust { c, delta } => Box::new(RandomDegreeSet::new(
                &SolitonDistribution::RobustSoliton {
                    dmax: params.num_active(),
                    c: *c,
                    delta: *delta,
                },
                params,
            )),
            DegreeSetConfig::RandomIdealTruncated { dmax } => Box::new(RandomDegreeSet::new(
                &SolitonDistribution::IdealSoliton {
                    dmax: *dmax.min(&params.num_active()),
                },
                params,
            )),
            DegreeSetConfig::TriangularIdealTruncated { dmax } => {
                Box::new(RandomTriangularDegreeSet::new(
                    &SolitonDistribution::IdealSoliton {
                        dmax: *dmax.min(&params.num_active()),
                    },
                    params,
                ))
            }
            DegreeSetConfig::ByName {
                name,
                degree_params,
            } => create_degree_set_by_name(name, degree_params, params),
        }
    }
}

/// Create degree set by name with parameters
pub fn create_degree_set_by_name(
    name: &str,
    degree_params: &[f64],
    code_params: &CodeParams,
) -> Box<dyn DegreeSetGenerator> {
    match name.to_lowercase().as_str() {
        "random_ideal" => Box::new(RandomDegreeSet::new(
            &SolitonDistribution::IdealSoliton {
                dmax: code_params.num_active(),
            },
            code_params,
        )),
        "random_robust" => {
            let c = degree_params.first().copied().unwrap_or(0.1);
            let delta = degree_params.get(1).copied().unwrap_or(0.5);
            Box::new(RandomDegreeSet::new(
                &SolitonDistribution::RobustSoliton {
                    dmax: code_params.num_active(),
                    c,
                    delta,
                },
                code_params,
            ))
        }
        "random_ideal_truncated" => {
            let dmax = degree_params.first().map(|&x| x as usize).unwrap_or(30);
            Box::new(RandomDegreeSet::new(
                &SolitonDistribution::IdealSoliton {
                    dmax: dmax.min(code_params.num_active()),
                },
                code_params,
            ))
        }
        "triangular_ideal_truncated" => {
            let dmax = degree_params.first().map(|&x| x as usize).unwrap_or(30);
            Box::new(RandomTriangularDegreeSet::new(
                &SolitonDistribution::IdealSoliton {
                    dmax: dmax.min(code_params.num_active()),
                },
                code_params,
            ))
        }
        "deterministic" => Box::new(DeterministicDegreeSet::new(code_params)),
        "fixed_degree" => {
            let degree = degree_params.first().map(|&x| x as usize).unwrap_or(3);
            Box::new(FixedDegreeSet::new(code_params, degree))
        }
        _ => panic!("Unknown degree set generator: {}", name),
    }
}

/*
/// Degree set type
#[derive(Clone, Debug)]
pub enum DegreeSetType {
    RandomIdeal,
    RandomRobust { c: f64, delta: f64 },
    RandomIdealTruncated { dmax: usize },
    TriangularIdealTruncated { dmax: usize },
}

impl DegreeSetType {
    pub fn new(degree_set_type: DegreeSetType, params: CodeParams) -> Box<dyn DegreeSetGenerator> {
        match degree_set_type {
            DegreeSetType::RandomIdeal => Box::new(
                RandomDegreeSet::new(
                    DegreeDistribution::IdealSoliton { dmax: params.num_active() },
                    params)
            ),
            DegreeSetType::RandomIdealTruncated { dmax } => Box::new(
                RandomDegreeSet::new(
                    DegreeDistribution::IdealSoliton { dmax: dmax.min(params.num_active()) },
                    params)
            ),
            DegreeSetType::RandomRobust { c, delta } => Box::new(
                RandomDegreeSet::new(
                    DegreeDistribution::RobustSoliton { dmax: params.num_active(), c, delta },
                    params)
            ),
            DegreeSetType::TriangularIdealTruncated { dmax } => Box::new(
                RandomTriangularDegreeSet::new(
                    DegreeDistribution::IdealSoliton { dmax: dmax.min(params.num_active()) },
                    params)
            )
        }
    }
}
*/
