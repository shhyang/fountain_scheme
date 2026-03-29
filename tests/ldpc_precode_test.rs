// Copyright (c) 2025 Shenghao Yang. All rights reserved.
// Licensed under the MIT License. See LICENSE-MIT for details.

use fountain_engine::decoder::Decoder;
use fountain_engine::encoder::Encoder;
use fountain_engine::traits::*;
use fountain_engine::types::*;
use fountain_scheme::validation::pseudo_rand::XorShift64;
use fountain_scheme::LDPCLTCode;
use fountain_utility::VecDataOperater;

/// LT encoding with R10-style LDPC precode: round-trip through encoder/decoder.
#[test]
fn test_ldpc_precode() {
    let k = 23;
    let t = 3;

    let mut message_vectors = vec![vec![0u8; t]; k];
    message_vectors[0][0] = 0;
    message_vectors[1][0] = 1;
    message_vectors[2][0] = 0;
    message_vectors[3][0] = 1;

    let mut vec_data_operater = VecDataOperater::new(t);
    message_vectors
        .iter()
        .enumerate()
        .for_each(|(i, vector)| vec_data_operater.insert_vector(vector, i));

    let config = LDPCLTCode::new_with_ideal_soliton(k, XorShift64::new(0x00C0_FFEE));
    let params = config.get_params();

    let mut encoder = Encoder::new_with_operator(config.clone(), Box::new(vec_data_operater));

    // Fountain symbols use coded_id >= num_total() (indices below that are message / LDPC / HDPC slots).
    let coded_ids = params.num_total()..params.num_total() + k * 3;
    for coded_id in coded_ids.clone() {
        encoder.encode_coded_vector(coded_id);
    }

    let mut decoder = Decoder::new_with_operator(config.clone(), Box::new(VecDataOperater::new(t)));

    let mut last = DecodeStatus::NotDecoded;
    for coded_id in coded_ids.clone() {
        let coded_vector = encoder.manager.get_coded_vector(coded_id);
        last = decoder.add_coded_vector(coded_id, &coded_vector);
        if last == DecodeStatus::Decoded {
            break;
        }
    }
    assert_eq!(
        last,
        DecodeStatus::Decoded,
        "LDPC-LT decode should finish within the emitted fountain symbols"
    );

    let dec_data_operater = decoder.manager.move_operator();
    for i in 0..k {
        assert_eq!(&dec_data_operater.get_vector(i), &message_vectors[i]);
    }
}
