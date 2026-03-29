// Copyright (c) 2025 Shenghao Yang. All rights reserved.
// Licensed under the MIT License. See LICENSE-MIT for details.


#[derive(Debug, Clone, Copy)]
pub enum SolitonDistribution {
    /// Ideal soliton distribution with maximum degree `dmax`.
    IdealSoliton {
        /// Maximum degree (truncation point).
        dmax: usize,
    },
    /// Robust soliton distribution parameterised by `c` (spike scale) and `delta` (failure bound).
    RobustSoliton {
        /// Maximum degree (truncation point).
        dmax: usize,
        /// Spike scale constant controlling the extra probability mass near degree `k/S`.
        c: f64,
        /// Target failure probability bound.
        delta: f64,
    },
}

impl SolitonDistribution {
    /// Computes and normalises the probability density function (indexed by degree).
    pub fn pdf(&self) -> Vec<f64> {
        let mut pdf = match self {
            SolitonDistribution::IdealSoliton { dmax } => self.ideal_soliton_pdf(*dmax),
            SolitonDistribution::RobustSoliton { dmax, c, delta } => {
                self.robust_soliton_pdf(*dmax, *c, *delta)
            }
        };
        let sum: f64 = pdf.iter().sum();
        if sum <= 0.0 {
            panic!("Invalid degree distribution: sum = {sum} is not positive");
        }
        pdf.iter_mut().for_each(|p| *p /= sum);
        pdf
    }

    fn ideal_soliton_pdf(&self, dmax: usize) -> Vec<f64> {
        let mut pdf = vec![0.0; dmax + 1];
        if dmax > 0 {
            pdf[1] = 1.0 / dmax as f64;
            for d in 2..=dmax {
                pdf[d] = 1.0 / (d as f64 * (d as f64 - 1.0));
            }
        }
        pdf
    }

    fn robust_soliton_pdf(&self, dmax: usize, c: f64, delta: f64) -> Vec<f64> {
        let ideal_pdf = self.ideal_soliton_pdf(dmax);
        if dmax == 0 {
            return ideal_pdf;
        }

        let kf = dmax as f64;
        let r = (kf.sqrt() * c * (kf / delta).ln()).max(1.0);
        let upper = (kf / r).ceil() as usize;

        let mut tau = vec![0.0; dmax + 1];
        for i in 1..upper.min(dmax) {
            tau[i] = r / (i as f64 * kf);
        }
        if upper <= dmax {
            tau[upper] = r * (r / delta).ln() / kf;
        }

        let beta: f64 = tau.iter().sum();
        let mut pdf = vec![0.0; dmax + 1];
        for i in 1..=dmax {
            pdf[i] = (ideal_pdf[i] + tau[i]) / (1.0 + beta);
        }
        pdf
    }
}
