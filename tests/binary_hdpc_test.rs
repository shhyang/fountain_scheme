// Copyright (c) 2025 Shenghao Yang. All rights reserved.
// Licensed under the MIT License. See LICENSE-MIT for details.

use fountain_engine::decoder::Decoder;
use fountain_engine::encoder::Encoder;
use fountain_engine::traits::*;
use fountain_engine::types::*;
use fountain_scheme::degree_sets::asymp_optimal_cdf;
use fountain_scheme::validation::pseudo_rand::XorShift64;
use fountain_scheme::BinaryHDPCLTCode;
use fountain_utility::VecDataOperater;

/// Ordinary mode via [`BinaryHDPCLTCode::new_with_ideal_soliton`] 
#[test]
fn test_binary_hdpc_lt_ordinary_from_ideal_soliton() {
    for k in 12..21 {
        for offset in 0..20 {
            let config = BinaryHDPCLTCode::new_with_ideal_soliton(k, XorShift64::new(0x00C0_FFEE));
            test_binary_hdpc_lt_roundtrip(&config, k, offset);
        }
    }
}

/// Ordinary mode via explicit CDF from [`BinaryHDPCLTCode::new`] 
#[test]
fn test_binary_hdpc_lt_ordinary_from_new_cdf() {
    for k in 12..21 {
        for offset in 0..20 {
            let cdf = asymp_optimal_cdf(k, 2, k);
            let config = BinaryHDPCLTCode::new(k, cdf, XorShift64::new(0xBEEF));
            test_binary_hdpc_lt_roundtrip(&config, k, offset);
        }
    }
}

fn test_binary_hdpc_lt_roundtrip(config: &BinaryHDPCLTCode<XorShift64>, k: usize, offset: usize) {
    let t = 1;

    let mut message_vectors = vec![vec![0u8; t]; k];
    message_vectors[k - 1][0] = 1;
    message_vectors[0][0] = 1;

    let mut vec_data_operater = VecDataOperater::new(t);
    message_vectors
        .iter()
        .enumerate()
        .for_each(|(i, vector)| vec_data_operater.insert_vector(vector, i));

    assert_eq!(config.code_type(), CodeType::Ordinary);

    let params = config.get_params();

    let mut encoder = Encoder::new_with_operator(config.clone(), Box::new(vec_data_operater));

    let total_num = params.num_total();
    let coded_ids = total_num + offset..total_num + 3 * k + offset;

    for coded_id in coded_ids.clone() {
        encoder.encode_coded_vector(coded_id);
    }

    let mut decoder = Decoder::new_with_operator(config.clone(), Box::new(VecDataOperater::new(t)));

    for coded_id in coded_ids {
        let coded_vector = encoder.manager.get_coded_vector(coded_id);
        let result = decoder.add_coded_vector(coded_id, &coded_vector);
        match result {
            DecodeStatus::Decoded => break,
            DecodeStatus::NotDecoded => {}
        }
    }

    for i in 0..k {
        assert_eq!(&decoder.manager.get_data_vector(i), &message_vectors[i]);
    }
}
