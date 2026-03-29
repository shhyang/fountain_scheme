// Copyright (c) 2025 Shenghao Yang. All rights reserved.
// Licensed under the MIT License. See LICENSE-MIT for details.

//! LDPC-LT Code Example (New Implementation)
//!
//! This example demonstrates LDPC-LT codes using the new LDPCLTCode struct.
//! This provides a cleaner, more direct implementation compared to using CodeConfig.

use fountain_scheme::validation::pseudo_rand::XorShift64;
use fountain_scheme::LDPCLTCode;
use fountain_engine::*;
use fountain_utility::VecDataOperater;

fn main() {
    println!("=== LDPC-LT Code Example (New Implementation) ===");

    // Configuration
    let k = 21; // number of source symbols
    let symbol_size = 4; // bytes per symbol

    println!("Parameters:");
    println!("  k (source symbols): {}", k);
    println!("  symbol_size: {} bytes", symbol_size);

    // Create test data with random-like patterns
    let mut message_vectors = vec![vec![0u8; symbol_size]; k];
    for i in 0..k {
        for j in 0..symbol_size {
            message_vectors[i][j] = ((i * 7 + j * 13) % 256) as u8; // Create distinctive patterns
        }
    }

    println!(
        "\nCreated {} message vectors of {} bytes each",
        k, symbol_size
    );
    println!("Sample data patterns:");
    for i in 0..std::cmp::min(5, k) {
        println!("  Symbol {}: {:?}", i, message_vectors[i]);
    }
    if k > 5 {
        println!("  ...");
    }

    // Setup encoder with data operator
    let mut vec_data_operater = VecDataOperater::new(symbol_size);
    for (i, vector) in message_vectors.iter().enumerate() {
        vec_data_operater.insert_vector(vector, i);
    }

    let ldpc_lt_code = LDPCLTCode::new_with_ideal_soliton(k, XorShift64::new(0xBEE));
    let params = ldpc_lt_code.get_params();

    println!("\nLDPC-LT Code Configuration:");
    println!("  Code type: {:?}", ldpc_lt_code.code_type());
    println!("  Parameters:");
    println!("    k: {}", params.k);
    println!("    a: {}", params.a);
    println!("    l: {}", params.l);
    println!("    b: {}", params.b);
    println!("    h: {}", params.h);
    println!("    Total symbols: {}", params.num_total());
    println!("    LDPC symbols: {}", params.l);

    // Create encoder
    let mut encoder = Encoder::new_with_operator(ldpc_lt_code.clone(), Box::new(vec_data_operater));

    // Generate coded symbols - LDPC-LT codes are more efficient than pure LT
    let num_coded_symbols = k * 2 + 20; // 2x overhead + buffer
    let start_id = params.num_total(); // Start after source + LDPC symbols
    let coded_ids: Vec<usize> = (start_id..start_id + num_coded_symbols).collect();

    let overhead_percent = ((num_coded_symbols as f64 / k as f64 - 1.0) * 100.0) as u32;
    println!(
        "\nGenerating {} coded symbols ({}% overhead)...",
        num_coded_symbols, overhead_percent
    );
    println!("  Starting coded symbol ID: {}", start_id);
    println!(
        "  Coded symbol range: {} to {}",
        start_id,
        start_id + num_coded_symbols - 1
    );

    for coded_id in &coded_ids {
        encoder.encode_coded_vector(*coded_id);
    }

    println!("Encoding complete!");

    // Simulate different loss scenarios
    let loss_scenarios = vec![0.0, 0.1, 0.2, 0.3]; // 0%, 10%, 20%, 30% loss

    for loss_rate in loss_scenarios {
        println!(
            "\n--- Testing {}% packet loss ---",
            (loss_rate * 100.0) as u32
        );

        // Simulate packet loss
        let num_received = ((num_coded_symbols as f64) * (1.0 - loss_rate)) as usize;
        let received_ids: Vec<usize> = coded_ids.iter().take(num_received).cloned().collect();

        println!(
            "Received {} out of {} coded symbols",
            num_received, num_coded_symbols
        );

        // Setup decoder
        let mut decoder = Decoder::new_with_operator(ldpc_lt_code.clone(), Box::new(VecDataOperater::new(symbol_size)));

        println!("Starting LDPC-LT decoding...");

        // Add received symbols to decoder
        let mut decoded_successfully = false;
        let mut symbols_used = 0;
        for (i, coded_id) in received_ids.iter().enumerate() {
            let coded_vector = encoder.manager.get_coded_vector(*coded_id);
            let status = decoder.add_coded_vector(*coded_id, &coded_vector);

            match status {
                DecodeStatus::Decoded => {
                    symbols_used = i + 1;
                    println!("✅ Successfully decoded after {} symbols!", symbols_used);
                    decoded_successfully = true;
                    break;
                }
                DecodeStatus::NotDecoded => {
                    if i % 3 == 0 && i > 0 {
                        // Print progress every 3 symbols
                        println!("  Progress: {}/{} symbols processed", i + 1, num_received);
                    }
                }
            }
        }

        if decoded_successfully {
            // Verify decoded data matches original
            let mut all_correct = true;
            for i in 0..k {
                let original = &message_vectors[i];
                let decoded = decoder.manager.get_data_vector(i);

                if *original != decoded {
                    all_correct = false;
                    break;
                }
            }

            if all_correct {
                let coded_overhead = ((symbols_used as f64 / k as f64 - 1.0) * 100.0) as u32;
                println!("✅ All symbols decoded correctly!");
                println!(
                    "✅ LDPC-LT code successfully handled {}% packet loss with {}% coded overhead",
                    (loss_rate * 100.0) as u32,
                    coded_overhead
                );
            } else {
                println!("❌ Decoding failed - some symbols are incorrect");
            }
        } else {
            println!(
                "❌ Decoding failed - insufficient symbols for {}% loss",
                (loss_rate * 100.0) as u32
            );
        }
    }

    println!("\n=== Example Complete ===");
}
