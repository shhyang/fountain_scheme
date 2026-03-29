// Copyright (c) 2025 Shenghao Yang. All rights reserved.
// Licensed under the MIT License. See LICENSE-MIT for details.

//! HDPC-LT Code Example
//!
//! This example demonstrates [`HDPCLTCode`]: R10-style LDPC + default cyclic HDPC precoding
//! and Raptor-style degree sets.

use fountain_scheme::validation::pseudo_rand::XorShift64;
use fountain_scheme::HDPCLTCode;
use fountain_utility::{
    test_code_scheme_multiple, test_code_scheme_with_data_vectors, TestResult, TestStatistics,
};

/// Compute standard deviation of decoding overhead (coded vectors used - k) for successful runs.
fn overhead_std_dev(k: usize, results: &[TestResult]) -> f64 {
    let overheads: Vec<f64> = results
        .iter()
        .filter(|r| r.num_mismatches == 0)
        .map(|r| r.decoding_metrics.num_coded_vectors as f64 - k as f64)
        .collect();
    if overheads.len() < 2 {
        return 0.0;
    }
    let mean = overheads.iter().sum::<f64>() / overheads.len() as f64;
    let variance = overheads
        .iter()
        .map(|x| (x - mean).powi(2))
        .sum::<f64>()
        / (overheads.len() - 1) as f64;
    variance.sqrt()
}

fn main() {
    println!("=== HDPC-LT Code Example ===\n");

    let k = 50;
    let data_vector_length = 4;
    let num_coded_vectors = 3 * k;

    println!("Parameters:");
    println!("  k (source symbols): {}", k);
    println!("  data_vector_length: {} bytes", data_vector_length);
    println!("  num_coded_symbols: {}", num_coded_vectors);

    println!("\n--- Testing HDPC-LT (R10 LDPC + default HDPC) ---");
    let code1 = HDPCLTCode::new_with_ideal_soliton(k, XorShift64::new(0xC0FFEE));
    let results1 = test_code_scheme_multiple(10, &code1, k, num_coded_vectors);
    let (num_runs1, successful_runs1, success_rate1) = TestStatistics::success_rate(&results1);
    let overhead_stats1 = TestStatistics::overhead_stats(k, &results1);
    let std_dev1 = overhead_std_dev(k, &results1);
    println!(
        "  Success: {}/{} ({:.1}%), mean overhead: {:.3}, std dev: {:.3}",
        successful_runs1, num_runs1, success_rate1 * 100.0, overhead_stats1.mean, std_dev1
    );

    // Test 2: HDPC-LT with message vectors for verification
    println!("\n--- Testing HDPC-LT with Message Vectors ---");
    let code2 = HDPCLTCode::new_with_ideal_soliton(k, XorShift64::new(0xBEEF));
    let mut results2 = Vec::new();
    for i in 0..5 {
        if i % 2 == 0 && i > 0 {
            println!("  Completed {} runs...", i);
        }
        let result = test_code_scheme_with_data_vectors(&code2, k, data_vector_length, num_coded_vectors);
        results2.push(result);
    }
    let (num_runs2, successful_runs2, success_rate2) = TestStatistics::success_rate(&results2);
    let overhead_stats2 = TestStatistics::overhead_stats(k, &results2);
    let std_dev2 = overhead_std_dev(k, &results2);
    println!(
        "  Success: {}/{} ({:.1}%), mean overhead: {:.3}, std dev: {:.3}",
        successful_runs2, num_runs2, success_rate2 * 100.0, overhead_stats2.mean, std_dev2
    );

    // Comparative analysis
    println!("\n=== Comparative Analysis ===");
    println!(
        "{:<40} {:>8} {:>8} {:>10} {:>10}",
        "Configuration", "Success%", "Runs", "Overhead", "Std Dev"
    );
    println!("{}", "-".repeat(80));
    println!(
        "{:<40} {:>7.1}% {:>8} {:>9.3} {:>9.3}",
        "HDPC-LT (R10+default HDPC)",
        success_rate1 * 100.0,
        num_runs1,
        overhead_stats1.mean,
        std_dev1
    );
    println!(
        "{:<40} {:>7.1}% {:>8} {:>9.3} {:>9.3}",
        "HDPC-LT (R10+default HDPC) w/ vectors",
        success_rate2 * 100.0,
        num_runs2,
        overhead_stats2.mean,
        std_dev2
    );

    // Best configuration by mean overhead (among those with at least one success)
    let configs = [
        ("HDPC-LT (R10+default HDPC)", success_rate1, num_runs1, &overhead_stats1, std_dev1),
        ("HDPC-LT (R10+default HDPC) w/ vectors", success_rate2, num_runs2, &overhead_stats2, std_dev2),
    ];
    let best = configs
        .iter()
        .filter(|(_, rate, _, _, _)| *rate > 0.0)
        .min_by(|(_, _, _, a, _), (_, _, _, b, _)| a.mean.partial_cmp(&b.mean).unwrap());

    if let Some((name, rate, _num_runs, stats, std_dev)) = best {
        println!("\n=== Best Configuration ===");
        println!("Configuration: {}", name);
        println!("  Success rate: {:.1}%", rate * 100.0);
        println!("  Mean overhead: {:.3}", stats.mean);
        println!("  Std deviation: {:.3}", std_dev);
        println!("  Min overhead: {:.3}", stats.min);
        println!("  Max overhead: {:.3}", stats.max);
    }

    println!("\n=== Example Complete ===");
}
