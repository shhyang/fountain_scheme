# Changelog

All notable changes to the published `fountain_scheme` crate are documented here.

## [1.3.0] - 2026-06-08

### Changed

- **`fountain_engine` dependency** bumped to **1.3.1** (native padding: `new_with_num_source`, `DataManager::set_num_source`, `Solver::install_padding`).
- **`fountain_utility` dependency** bumped to **1.3.0** (`padding_codec` wrappers for scheme-agnostic padding).
- Integration tests and examples updated for `Encoder` / `Decoder` APIs that take `&CodeScheme` references.
- **RaptorQ HDPC** (`precodes/hdpc.rs`): `mul_binary` uses row-buffered `rq_mul_binary_from_rows` for fewer allocations during inactivation GE.
- Scheme sources synced from monorepo (`degree_sets`, `parameters`, `validation`, LT/LDPC/HDPC scheme types).

### Notes

- Published against the **legacy solver** `fountain_engine` surface (no `BinaryMatrix` / `mul_binary_packed` fast path). Binary HDPC uses dense `mul_binary` only; behavior matches prior releases.
