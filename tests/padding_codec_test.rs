// Copyright (c) 2025 Shenghao Yang. All rights reserved.
// Licensed under the MIT License. See LICENSE-MIT for details.

//! Integration tests for [`fountain_utility::padding_codec`] with a real scheme.

use fountain_engine::traits::{CodeScheme, DataOperator};
use fountain_engine::types::{CodeType, DecodeStatus};
use fountain_scheme::HDPCLTCode;
use fountain_scheme::validation::pseudo_rand::XorShift64;
use fountain_utility::VecDataOperater;
use fountain_utility::padding_codec::{BlockSizePolicy, PaddedDecoder, PaddedEncoder, WithSourceK};

#[test]
fn ordinary_padded_round_trip_with_source_k() {
    let block_k = 12;
    let source_k = 10;
    let t = 8;
    let inner = HDPCLTCode::new_with_ideal_soliton(block_k, XorShift64::new(42));
    let scheme = WithSourceK::new(inner, source_k);
    assert_eq!(scheme.code_type(), CodeType::Ordinary);
    assert_eq!(scheme.num_padding(), 2);

    let mut messages = vec![vec![0u8; t]; source_k];
    for (i, row) in messages.iter_mut().enumerate() {
        for (j, byte) in row.iter_mut().enumerate() {
            *byte = ((i + 1) * (j + 3) % 251) as u8;
        }
    }

    let mut enc_op = VecDataOperater::new(t);
    for (i, v) in messages.iter().enumerate() {
        enc_op.insert_vector(v, i);
    }

    let mut encoder = PaddedEncoder::new_with_operator(scheme.clone(), Box::new(enc_op), t);
    let params = scheme.get_params();
    let num_coded = source_k + source_k / 2 + scheme.num_padding();
    let coded_ids: Vec<usize> = (params.num_total()..params.num_total() + num_coded).collect();

    let mut coded_payload = Vec::new();
    for &coded_id in &coded_ids {
        if let Some(data_id) = encoder.encode_coded_vector(coded_id) {
            coded_payload.push((coded_id, encoder.get_data_vector(data_id).to_vec()));
        }
    }

    let dec_op = VecDataOperater::new(t);
    let mut decoder = PaddedDecoder::new_with_operator(scheme, Box::new(dec_op), t);

    let mut decoded = false;
    for (coded_id, payload) in coded_payload {
        if decoder.add_coded_vector(coded_id, &payload) == DecodeStatus::Decoded {
            decoded = true;
            break;
        }
    }
    assert!(decoded, "ordinary padded decode should succeed");

    let dec_op = decoder.inner.manager.move_operator();
    for (i, expected) in messages.iter().enumerate() {
        assert_eq!(dec_op.get_vector(i), expected.as_slice(), "message {i}");
    }
}
