// Copyright (c) 2025 Shenghao Yang. All rights reserved.
// Licensed under the MIT License. See LICENSE-MIT for details.

use fountain_engine::decoder::Decoder;
use fountain_engine::encoder::Encoder;
use fountain_engine::traits::*;
use fountain_engine::types::*;
use fountain_scheme::validation::pseudo_rand::XorShift64;
use fountain_scheme::HDPCLTCode;
use fountain_utility::VecDataOperater;

/// [`HDPCLTCode::as_systematic`] must flip [`CodeScheme::code_type`] (regression guard).
#[test]
fn test_hdpc_lt_as_systematic_sets_code_type() {
    let k = 12;
    let mut code = HDPCLTCode::new_with_ideal_soliton(k, XorShift64::new(1));
    assert_eq!(code.code_type(), CodeType::Ordinary);
    code.as_systematic();
    assert_eq!(code.code_type(), CodeType::Systematic);
}

/// Full systematic encode/decode with LDPC + HDPC precode (RQ-style HDPC / RQLDPC).
///
/// Same flow as `tests/systematic_ldpc_test.rs` (`test_systematic_ldpc_lt_roundtrip`):
/// systematic indices `0..k`, then extra fountain symbols from `num_total()`.
#[test]
fn test_systematic_hdpc_lt_roundtrip() {
    let k = 36;
    let t = 3;

    let mut message_vectors = vec![vec![0u8; t]; k];
    message_vectors[0][0] = 1;

    let mut vec_data_operater = VecDataOperater::new(t);
    message_vectors
        .iter()
        .enumerate()
        .for_each(|(i, vector)| vec_data_operater.insert_vector(vector, i));

    let mut config = HDPCLTCode::new_with_ideal_soliton(k, XorShift64::new(0x00C0_FFEE));
    config.as_systematic();

    let params = config.get_params();

    let mut encoder = Encoder::new_with_operator(config.clone(), Box::new(vec_data_operater));

    for coded_id in 0..k {
        encoder.encode_coded_vector(coded_id);
    }
    let fountain_ids = params.num_total()..params.num_total() + k * 3;
    for coded_id in fountain_ids.clone() {
        encoder.encode_coded_vector(coded_id);
    }

    let mut decoder = Decoder::new_with_operator(config.clone(), Box::new(VecDataOperater::new(t)));

    let mut last = DecodeStatus::NotDecoded;
    for coded_id in 0..k {
        let v = encoder.manager.get_coded_vector(coded_id);
        last = decoder.add_coded_vector(coded_id, &v);
        if last == DecodeStatus::Decoded {
            break;
        }
    }
    if last != DecodeStatus::Decoded {
        for coded_id in fountain_ids {
            let v = encoder.manager.get_coded_vector(coded_id);
            last = decoder.add_coded_vector(coded_id, &v);
            if last == DecodeStatus::Decoded {
                break;
            }
        }
    }

    assert_eq!(last, DecodeStatus::Decoded);

    let dec = decoder.manager.move_operator();
    for i in 0..k {
        assert_eq!(&dec.get_vector(i), &message_vectors[i]);
    }
}
