# C147 saved-result independent audit

Verdict: PASS for the seven saved exact density stages of the two frozen source cases. This is a saved-output audit, not a replay of the author generator and not a kernel/moment/C6 computation.

The reviewed source result is `/Users/zhengzhihao/Desktop/T/Math-Research-cube6-sol-localmax-audit/cube-dim-6-bounds/upper/round147/pair_pushforward/next_result.json`, SHA256 `1233ffbd6c9b093f3e99715affe3ef274e4d7fff374e23f9f3aa4c2c9889acd8`, 9941 bytes. Its exit receipt SHA256 is `c5d44067b73d8cddfac3420aa6290b0406db7b257443cfe2ec5456a798b5063f`. The source result, exit receipt, start record and exact-hash release agree. The main release carries the prior independent review manifest and the external literal endpoint budget34.

## Actual independent batch

`registration.json` was written before importing SymPy or doing symbolic reconstruction. One batch ran with SymPy1.14.0, nice10 and a one-thread environment. It read the saved coefficients through a restricted arithmetic parser, not eval/import of the author program. It reconstructed all seven density-sector Laplace integrals directly from their saved rational rays, covectors and polynomial coefficients.

For every stage the following passed independently:

- The case, zero-column count and unchanged six source columns match the original pinned C130 cone table and the frozen two-case packet. Active prefixes are independently reconstructed from the first independent pivot and remaining original columns.
- The saved cells are exactly the consecutive distinct source directions, correctly oriented and normalized, covering the prefix column cone with disjoint interiors and only measure-zero boundary overlap.
- All coefficients are rational; all nonzero polynomial degrees obey the stage bound m-2; equal covectors are not duplicated within a cell.
- Each saved exponential covector has strictly positive dot product with both of its support boundary rays. Consequently every saved signed polynomial-exponential term is absolutely integrable on its own support, including at Laplace parameters0. Negative coordinates of a covector cause no defect on its designated cone.
- The reconstructed rational Laplace transform minus the original independent-exponential product has exactly zero numerator. The separately parsed saved product also agrees with the source product. The reconstructed total mass is exactly1.

| Case | Stage | Sectors | Terms | Maximum degree | Minimum boundary-ray decay |
|---|---:|---:|---:|---:|---:|
| m3 body(9,16) | 2 | 1 | 1 | 0 | 30 |
| m3 body(9,16) | 3 | 2 | 4 | 0 | 30 |
| m6 body(0,1) | 2 | 1 | 1 | 0 | 10 |
| m6 body(0,1) | 3 | 1 | 2 | 0 | 10 |
| m6 body(0,1) | 4 | 1 | 2 | 1 | 10 |
| m6 body(0,1) | 5 | 1 | 2 | 2 | 10 |
| m6 body(0,1) | 6 | 1 | 3 | 2 | 10 |

The final m3 density covers the nonnegative quadrant in two sectors separated by A=B. The final m6 density is supported on B<=A<=2B. The complete saved signed density has the same exact Laplace transform as the positive input-exponential pushforward law. Together with absolute integrability and uniqueness of the finite-measure Laplace transform, this identifies the saved density with that probability law; the earlier independently audited convolution proof supplies the same positivity mechanism. No signed summand is assumed positive.

Source pins were checked again after the batch. Independent reconstruction wall0.576681792 seconds, CPU0.494171 seconds. No failures or reruns occurred.

## Bookkeeping and scope limits

The source program reports16 under `raw_endpoint_contributions`. Its resonance branch counts one endpoint difference instead of two literal endpoints. This audit does not call16 an actual physical endpoint count and does not claim an exact physical count. The prior independent frozen-case budget remains valid: at most4 literal endpoints for m3, at most2+4+8+16=30 for m6, at most34 combined. This is below the released cap100000. The convention is already disclosed by the main release; it does not affect any exact density identity.

The result/exit flags report zero kernel/native/grid calls and no new moment or C6 bound. Source/gate inspection and this saved audit are consistent with those flags. This auditor performed zero new density generations, zero old m2 toy repetitions, zero kernel or special-function value evaluations, zero native queries and zero old6D replays. The author generator was neither imported nor executed. No149 patterns or other new scientific cases were generated. All writes belong only to this clone's `cube-dim-6-bounds/audit/upper147_saved`; nothing was staged, committed or pushed.
