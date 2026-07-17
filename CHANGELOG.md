# Changelog

All notable changes to the published `fountain_scheme` crate are documented here.

## [2.0.0] - 2026-07-17

### Changed

- **`fountain_engine` dependency** bumped to **2.0.0** (SystemSolver-based engine).
- **`fountain_utility` dependency** bumped to **2.0.0**.
- **HDPC precodes:** `mul_data` paths use `add_one_to_vector` / `broadcast_add_owned`; validation uses `ensure_zero_one`.
- **Binary HDPC:** adds `mul_binary_packed` for the engine v2 packed GF(2) inactive path (in addition to dense `mul_binary`).

## [1.3.0] - 2026-06-08

### Changed

- **`fountain_engine` dependency** bumped to **1.3.1** (native padding: `new_with_num_source`, `DataManager::set_num_source`, `Solver::install_padding`).
- **`fountain_utility` dependency** bumped to **1.3.0** (`padding_codec` wrappers for scheme-agnostic padding).
- Integration tests and examples updated for `Encoder` / `Decoder` APIs that take `&CodeScheme` references.
- **RaptorQ HDPC** (`precodes/hdpc.rs`): `mul_binary` uses row-buffered `rq_mul_binary_from_rows` for fewer allocations during inactivation GE.
- Scheme sources synced from monorepo (`degree_sets`, `parameters`, `validation`, LT/LDPC/HDPC scheme types).

### Notes

- Published against the then-current `fountain_engine` 1.x surface; binary HDPC used dense `mul_binary` only.
