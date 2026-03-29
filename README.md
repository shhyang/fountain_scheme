# fountain_scheme

Configurable fountain code schemes built on `fountain_engine`.

## Included Schemes

- **LT Code**: Pure Luby Transform codes with ideal soliton and robust soliton degree distributions.
- **LDPC-LT Code**: LT codes with LDPC precode for improved decoding performance.
- **HDPC-LT Code**: LT codes with HDPC (High-Density Parity Check) precode.
- **Systematic Codes**: Systematic variants of LDPC-LT and HDPC codes, where source symbols appear unmodified in the output.

## Quick Start

Use the code scheme types (each implements `fountain_engine::traits::CodeScheme`). Pass any of them to `fountain_engine::Encoder` and `fountain_engine::Decoder`.

```rust
use fountain_scheme::validation::pseudo_rand::XorShift64;
use fountain_scheme::{BinaryHDPCLTCode, HDPCLTCode, LDPCLTCode, RandomLTCode};

let k = 100usize;
let rng = XorShift64::new(1);

let _lt = RandomLTCode::new_from_robust_soliton(k, rng.clone());
let _ldpc_lt = LDPCLTCode::new_with_ideal_soliton(k, rng.clone());
let _hdpc_lt = HDPCLTCode::new_with_ideal_soliton(k, rng.clone());
let mut sys = HDPCLTCode::new_with_ideal_soliton(k, rng.clone());
sys.as_systematic();
let _binary_hdpc = BinaryHDPCLTCode::new_with_ideal_soliton(k, rng);
```

## License

MIT License. See [LICENSE-MIT](LICENSE-MIT).

Copyright (c) 2025 Shenghao Yang. All rights reserved.
