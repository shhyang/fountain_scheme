// Copyright (c) 2025 Shenghao Yang. All rights reserved.
// Licensed under the MIT License. See LICENSE-MIT for details.

use fountain_engine::decoder::Decoder;
use fountain_engine::encoder::Encoder;
use fountain_engine::traits::*;
use fountain_engine::types::*;
use fountain_scheme::validation::pseudo_rand::XorShift64;
use fountain_scheme::{RandomLTCode, DeterministicLTCode};
use fountain_utility::VecDataOperater;

/// We test the library for LT code with no precode.
#[test]
fn test_lt_ideal() {
    let k = 15; // number of message vectors
    let t = 3; // vector length

    let mut message_vectors = vec![vec![0u8; t]; k];
    message_vectors[0][0] = 1;

    // Create VecDataOperater for actual data storage
    let mut vec_data_operater = VecDataOperater::new(t);
    message_vectors
        .iter()
        .enumerate()
        .for_each(|(i, vector)| vec_data_operater.insert_vector(vector, i));

    // encoder
    let config = RandomLTCode::new_from_robust_soliton(k, XorShift64::new(0x00C0_FFEE));
    // Record encoding operations in basic_manager
    let mut encoder = Encoder::new_with_operator(config.clone(), Box::new(vec_data_operater));

    let coded_ids = k..k * 3;
    // encode the coded vectors
    for coded_id in coded_ids.clone() {
        encoder.encode_coded_vector(coded_id);
    }

    // decoding
    let mut decoder = Decoder::new_with_operator(config.clone(), Box::new(VecDataOperater::new(t)));

    for coded_id in coded_ids.clone() {
        let coded_vector = encoder.manager.get_coded_vector(coded_id);
        let result = decoder.add_coded_vector(coded_id, &coded_vector);
        match result {
            DecodeStatus::Decoded => {
                println!("coded_id: {coded_id}, status: {result:?}");
                break;
            }
            DecodeStatus::NotDecoded => {}
        }
    }
    let dec_data_operater = decoder.manager.move_operator();

    for i in 0..k {
        assert_eq!(&dec_data_operater.get_vector(i), &message_vectors[i]);
    }
}


/// We test the library for LT code with no precode.
#[test]
fn test_deterministic_lt_ideal() {
    let k = 15; // number of message vectors
    let t = 3; // vector length

    let mut message_vectors = vec![vec![0u8; t]; k];
    message_vectors[0][0] = 1;

    // Create VecDataOperater for actual data storage
    let mut vec_data_operater = VecDataOperater::new(t);
    message_vectors
        .iter()
        .enumerate()
        .for_each(|(i, vector)| vec_data_operater.insert_vector(vector, i));

    // encoder
    let config = DeterministicLTCode::new_from_robust_soliton(k);
    // Record encoding operations in basic_manager
    let mut encoder = Encoder::new_with_operator(config.clone(), Box::new(vec_data_operater));

    let coded_ids = k..k * 3;
    // encode the coded vectors
    for coded_id in coded_ids.clone() {
        encoder.encode_coded_vector(coded_id);
    }

    // decoding
    let mut decoder = Decoder::new_with_operator(config.clone(), Box::new(VecDataOperater::new(t)));

    for coded_id in coded_ids.clone() {
        let coded_vector = encoder.manager.get_coded_vector(coded_id);
        let result = decoder.add_coded_vector(coded_id, &coded_vector);
        match result {
            DecodeStatus::Decoded => {
                println!("coded_id: {coded_id}, status: {result:?}");
                break;
            }
            DecodeStatus::NotDecoded => {}
        }
    }
    let dec_data_operater = decoder.manager.move_operator();

    for i in 0..k {
        assert_eq!(&dec_data_operater.get_vector(i), &message_vectors[i]);
    }
}
