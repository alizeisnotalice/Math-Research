# R062 compute-facts review

**Finite-run facts: ACCEPT.** Read-only provenance and execution review of the two completed frontier batches and their separately generated exact certificates. No grid or integral was rerun, no certificate proof was rederived, and no author/main-clone files were changed.

## Batch 1: fixed pilot

- Registration `frontier_registration.json` SHA-256 `08ccea426c7663cab0b1140fe848da65ac7908ef4f03a5e47c19219dca238fb1` binds `frontier.py` SHA-256 `6ded8207232f07d0ea5364d5c03d3eb0b2b7ebe99f18326d4e98e844c0106c2a`; current source matches. Before starting workers, the source asserts that its hash equals the registered hash.
- `frontier.log` records start PID 91481, 128 tasks, 15 workers, and the matching registration hash, followed by completion at 0.784092 s with 4,096 tuples, 245,760 band evaluations, zero errors, and 3,098 eligible entries. The summary binds the same registration hash.
- All 128 chunk files exist. Each has a completed status, and its SHA-256 matches the corresponding summary receipt. The 128 unique `(i,J)` keys cover 32 `i` indices × four layer sizes `{16,32,64,128}`. The raw chunk schema has no source/code hash field; the registration pins source identity and the summary pins each chunk hash.
- The summary `top64` exactly matches merging the saved `top32` from every chunk and sorting the retained values. This validates the finite result merge, not the ideal model or a global extremum.
- Per-chunk CPU receipts sum to 7.853184 s across 15 distinct worker PIDs. The only monitor sample is at 0.54 s: 1051.2% aggregate CPU (about 10.51 cores), and 454,320 KiB combined RSS. This single snapshot is not a sustained-load measurement.

## Batch 2: frontier v2

- Registration `frontier_v2_registration.json` SHA-256 `9252ff3e77c6d3c7ad6ea4701e229bfd2b2033254001fd9fbababaff383f59be` binds `frontier_v2.py` SHA-256 `aa9da32b4f6d195bdc44f7e57013c6c5f8737aaf98003b5482119f281abec4ed`; current source matches and the run source guards against a mismatch before launching. Registration records 256 chunks and a zero old-pilot Fraction overlap.
- `frontier_v2.log` records start PID 91776, 256 tasks, 15 workers and the matching registration hash. Its completion event records 136.016310 s, 131,072 tuples, 50,331,648 band evaluations, zero errors and 84,808 eligible entries. The summary binds the same registration.
- All 256 chunk files exist, are completed, and each chunk hash matches its summary receipt. The 256 unique chunk indices cover `i=0,…,255`; each chunk contains both registered layer sizes 256 and 512. Raw chunk files do not include a source hash. `frontier_v2_summary.json`'s retained `top64` exactly matches merging and sorting the saved `top32` per chunk.
- The separate `frontier_v2_postrun_manifest.json` SHA-256 `1fbf1136884c8f9f7e4813f98abdc97a493b39a4bd500c03d45c4450812c90c2` associates the original registration, summary, current source hash and all 256 raw chunk hashes after the run. I checked every manifest chunk digest against the raw file and summary map; all match. It reports unchanged chunk files. This is explicitly post-run source association; it must not be described as a source hash embedded in each original raw receipt.
- Per-chunk CPU receipts sum to 1,846.221808 s across 15 distinct worker PIDs, about 13.57 core-equivalents averaged over the 136.02 s wall interval. The 27 `frontier_v2_monitor.jsonl` samples peak at 1458.8% aggregate CPU (about 14.588 cores) and 458,416 KiB combined RSS (about 448 MiB). The CPU peak is sampled, not sustained peak utilization; no evidence indicates memory exhaustion.

## Separate exact certificates and scope

- `candidate_certificate.json` SHA-256 `b889e35ef909bd979937752c67410bf3913eb76cdf7dc6c720556e093d4c0711` contains 128 fixed-layer rows (`J=128`, fixed `d=7/25`, and the listed fixed scalar parameters). Its `code_sha256` matches `candidate_cert.py` (`b51360a66b965ff732c579cbcf347068e73285174210f27f9f719d6144b86884`). It labels its scope “conditional scalar model; no operator theorem.” This review checks identity and scope only, not its arithmetic proof.
- `candidate_gate.json` SHA-256 `393aae74e54a2deae907e83b2b710786652ec0624fbc046bcf519be016caae9a` is a later post-run scalar gate. Its input hash matches the candidate certificate, and its source hash matches `candidate_gate.py` (`451c160f543283f03cfe44d2cbf16e02c0c51d9c54f5db3ff402a97a8da51d5d`). It records `gamma_threshold=429/2000`, `gamma_gt_gate=true`, and an exact positive margin. It is an add-on check and does not revise or retroactively alter the original candidate certificate.
- `frontier_v2_top_certificate.json` SHA-256 `ba081481b75d2c9eac25709f069fcf631d1ce0b2b9ef1bcdd8724f5b1411efd3` binds the v2 summary hash and its code hash matches `frontier_v2_top_cert.py` (`579c82a82faea04654a95c4e38b9746617397cb4767f761f03ed30301c0ce7c5`). It contains 512 layer rows for the retained top parameter `(i=183,k=125,J=512)` and reports positive/eligibility certificates for that finite conditional-model candidate. It does not establish global optimality over all parameters, and it is not a `C_n` theorem.
- Both frontier registrations explicitly scope the batches to a finite layered ideal model. The retained top entries, candidate certificate, gate, and v2 top certificate do not prove an operator bound or replace the still-needed R061 direct-sum proof in the mathematical argument.

## Hashes of primary receipts

`frontier_summary.json`: `8303a0aad22e13495d30354789e0e3a04d7dbc0f535072817f42ec9a5bd7e964`; `frontier_monitor.jsonl`: `01f5fd6d45b3ca0c8b8440a9b99e2dbd8a20666483b6230d62080165a581fd4d`.

`frontier_v2_summary.json`: `f33ae4461e784598329095ad4735806a5559794b937e1094dadc0417f2d767c2`; `frontier_v2_monitor.jsonl`: `1041cc8d18c2a92abe689578bc13297db9bf8dd304083b5b2d154e85cb78a998`; `frontier_v2_postrun_manifest.json`: `1fbf1136884c8f9f7e4813f98abdc97a493b39a4bd500c03d45c4450812c90c2`.
