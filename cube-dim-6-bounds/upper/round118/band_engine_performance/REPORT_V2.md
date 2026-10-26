# Round 118 precision calibration v2

## Outcome

The main-approved, task-owned precision-only C016 copy compiled at 128 bits and ran four registered saved point-window inputs sequentially. All four native processes exited 0 with empty stderr. For each output, the exact dyadic interval check confirms containment of the saved 512-bit interval for the band, all six gradients, `p_a`, `p_T`, and Euler. No 256-bit fallback was needed. This is finite calibration evidence for these four inputs, not a global certificate, theorem, or C6 bound. Production round117 remains on 512 bits.

The fourth v1 entry was replaced in v2 by the sealed round115 pilot1 lowtime input `B_[1/64,1/32]`, center `(1,7/8,7/8,7/8,1/4,1/4)`. Its saved 512-bit native SHA-256 is `e6ce0a7a87bca388e18c074d30d31e53eace990d8b26ca07582e51acec259a20`; receipt reports exit 0, empty stderr, and the full native output contains band, six gradients, `p_a`, `p_T`, and Euler. The v1 registration remains preserved and unrun.

## Per-case measurements

Times are measurements from these individual processes. Existing 512-bit `kernel_seconds` came from earlier, separate runs; no fair speedup ratio is claimed. Their old receipts do not provide comparable outer full-wall timing.

| Case | Window | 128-bit kernel seconds | Native process wall seconds | Child CPU seconds | Exit | Saved-512 containment |
|---|---|---:|---:|---:|---:|---|
| `C114.r0110101.a0.q009` | `[1/32,1/2]` | 2.90944 | 3.4303 | 2.9125 | 0 | pass, all 10 intervals |
| `C114.r011111.a0.q007` | `[1/32,1/2]` | 4.24042 | 4.2455 | 4.2375 | 0 | pass, all 10 intervals |
| `C114.r01011.a0.q003` | `[1/32,1/2]` | 3.28247 | 3.2871 | 3.2813 | 0 | pass, all 10 intervals |
| `C115.pilot1` | `[1/64,1/32]` | 3.42145 | 3.4261 | 3.4199 | 0 | pass, all 10 intervals |

“10 intervals” means band + six gradient coordinates + `p_a` + `p_T` + Euler. Comparison uses exact rational endpoints decoded from each dyadic mantissa/exponent pair. Per-field results and exact process records are in `run_receipt_v2.json`; complete stdout JSON outputs are under `runs_v2/`.

## Build, scope, and integrity

- The only source change is the registered compile-time `BAND_ENGINE_PRECISION_BITS` parameter; original C016 source, binary, fan input, FLINT/GMP/MPFR libraries, and header archive pins matched their recorded hashes before compilation and remain unchanged on the post-run check.
- 128-bit executable SHA-256: `d50e5822b3e3d6d40550386b0aaaf118569acd55d1214a3926f51cadcba3f9d5`. Compile exited 0; post-link path repair exited 0. Build timings are recorded in `build_receipt_128.json`.
- Four new native calls total; no 512-bit rerun, 256-bit compile/call, extra case, MATLAB, stage, commit, or push.
- `started_v2.json` is the exclusive-create pre-run record and was not modified. `approval_snapshot_v2.txt`, the finite `test_registration_v2.json`, saved references, source hashes, run receipt, and output hashes are recorded in `manifest_v2.json` / `final_receipt_v2.json`.
