// Copyright (c) 2025 Shenghao Yang. All rights reserved.
// Licensed under the MIT License. See LICENSE-MIT for details.

use fountain_engine::decoder::Decoder;
use fountain_engine::encoder::Encoder;
use fountain_engine::traits::*;
use fountain_engine::types::*;
use fountain_scheme::validation::pseudo_rand::XorShift64;
use fountain_scheme::HDPCLTCode;
use fountain_utility::VecDataOperater;

/// Integration test: ordinary HDPC-LT encode/decode with R10 precoding.
#[test]
fn test_hdpc_precode_1() {
    //test_hdpc_precode(259, 11);
    //test_hdpc_precode(21, 0);
    for k in 12..21 {
        for offset in 0..20 {
            test_hdpc_precode(k, offset);
        }
    }
}

fn test_hdpc_precode(k: usize, offset: usize) {
    let t = 1; // vector length

    //let message_vectors = gen_random_matrix(k, t, 8);
    let mut message_vectors = vec![vec![0u8; t]; k];
    message_vectors[k - 1][0] = 1;
    message_vectors[0][0] = 1;

    // Create VecDataOperater for actual data storage
    let mut vec_data_operater = VecDataOperater::new(t);
    message_vectors
        .iter()
        .enumerate()
        .for_each(|(i, vector)| vec_data_operater.insert_vector(vector, i));

    let config = HDPCLTCode::new_with_ideal_soliton(k, XorShift64::new(0x00C0_FFEE));
    let params = config.get_params();
    // println!("test of k: {}, offset: {}", k, offset);
    // println!("ldpc l: {}", params.l);
    // println!("b: {}", params.b);
    // println!("hdpc h: {}", params.h);
    // println!("num_pre_inactive: {}", params.num_pre_inactive());

    // Record encoding operations in basic_manager
    let mut encoder = Encoder::new_with_operator(config.clone(), Box::new(vec_data_operater));

    let total_num = params.num_total();
    let coded_ids = total_num + offset..total_num + 3 * k + offset;

    for coded_id in coded_ids.clone() {
        encoder.encode_coded_vector(coded_id);
    }

    // decoder
    let mut decoder = Decoder::new_with_operator(config.clone(), Box::new(VecDataOperater::new(t)));

    for coded_id in coded_ids {
        let coded_vector = encoder.manager.get_coded_vector(coded_id);
        let result = decoder.add_coded_vector(coded_id, &coded_vector);
        match result {
            DecodeStatus::Decoded => {
                //println!("decoded with coded_id: {}", coded_id);
                break;
            }
            DecodeStatus::NotDecoded => {}
        }
    }

    for i in 0..k {
        //println!("get_data_vector({}): {:?}", i, decoder.manager.get_data_vector(i));
        //println!("message_vectors[{}]: {:?}", i, message_vectors[i]);
        assert_eq!(&decoder.manager.get_data_vector(i), &message_vectors[i]);
    }
}
