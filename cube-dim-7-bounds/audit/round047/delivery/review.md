# R047 delivery fact review

**Disposition: FACTS ACCEPT.** Read-only review of producer evidence and archive identity; no quadrature was rerun, and no mathematical theorem audit was attempted. This review does not modify author or main-clone files.

## Registered run and retained data

- Producer source/report: `/Users/zhengzhihao/Desktop/T/cube-r030-compute/cube-dim-7-bounds/upper/round047/compute/`.
- `registration.json` SHA-256 `9ea30edad899d99cde9f6d409e23b8f2516ce1ecd580661622327dd5f638b4ce` pins `cross_face.py` SHA-256 `cd36c22d1f9dd78fd15a0f22e438dd7bb9f2bc7f8cf1fe6257dbad17a29a9551`, 15 workers, 80 dps, 4,096 triples per task, and the task stop condition. `summary.json` references the same registration hash and records 15 completed tasks. Its SHA-256 is `ba2057522d57e5e369bc89b07d650ed92e68104e1acdb9ec48ee5901d791be01`.
- The 15 original `case_*.json` raw files remain in this compute clone. Their per-task headers record 4,096 completed triples, 2,176 unique cross integrals, and 10,112 quadrature parts each. Cross-checking against archive-manifest raw hashes finds all 15 current raw files match their listed hashes.
- Totals reconcile: 15×4,096 = **61,440 ordered triples**; 15×2,176 = **32,640 unique cross-radius integrals**; 15×10,112 = **152,064 quadrature segments**. Aggregate wall is 61.0558 s, summed task CPU is 385.0261 s (about 6.31 core-equivalents averaged across elapsed wall time), with 15 distinct worker PIDs. Maximum per-worker RSS is 104,431,616 bytes (about 99.6 MiB).

## Numerical findings and scope

- The archived summary selects the min/max row from each of the 15 frozen raw cases. The lowest stored observation is `0.99999998901652094561531400997836` at `case_m32_d65536`, `u/d=1/4`, `theta1=theta2=-1`. The largest archived observation is the rounded value `1.0` in `case_m128_d256`, `u/d=1/4`, `theta1=theta2=1`. Both rows occur in the corresponding retained raw file. The producer notes that a displayed `1.0` at high `m` can reflect finite precision.
- The calculation uses 80-digit mpmath and forms `mu−1` directly, then uses `log1p` before multiplying by `n`; this is the stated numerical stabilization method. It avoids constructing an `n`-length vector. These implementation choices and the finite data support reproducibility of the diagnostic but do not prove a global inequality.
- `check_moments.json` records two supplementary 250-dps direct moment recomputations, at `m=32` and `m=512`, each agreeing with its stored row within about `4×10^-33`. This is a post-run numerical cross-check by the project producer, not an independent second-author audit, and not a rigorous interval enclosure.
- The registration stop was 4,096 triples per task (or stop on a numerical violation/invalid geometry). The sampled maxima and archive verifier cannot establish maximality over unrecorded triples beyond the registered grid. “No sampled value exceeded 4,096” is not a global bound and does not establish a `C_n` result or any new operator estimate.

## Resource samples

`monitor.jsonl` SHA-256 is `8cc42ffae01e766c11366c5f54e27667693e9febba32347f85cd9bd7d3035f1a`; it contains 13 samples. The peak CPU sample at 5.152 s is 1461.3% aggregate across 15 workers (about 14.61 cores) and RSS 555,200 KiB. CPU then falls: 1346.3% at 10.20 s; about 593–599% from 15.24–30.34 s; 585.6% at 35.38 s; 389.7% at 40.43 s; about 296% at 45–50.51 s; and 117.8% at 55.54 s, with the last 60.58 s sample at 100.0%. Monitored aggregate RSS peaks at 934,752 KiB (about 913 MiB). These samples show a pronounced declining-load tail; the brief 14.61-core peak is not 61 seconds of full-core loading, nor evidence that memory was exhausted.

## Original run versus post-run verification/archive

- The original-run identity is the registration above: its code hash is the `cross_face.py` hash and its per-task output hashes are in `summary.json` and raw case files.
- `check_moments.py`/`check_moments.json` are later supplemental verification artifacts, not part of the registered original run. Their script was timestamped after the original summary.
- `make_archive.py`, `archive_manifest.json`, `archive_summary.json.gz`, and `verify_archive.py`/receipt were created after the numerical run. The manifest source hash `357005990a7199f3370aea8dbeb9aacfb25835d9d23561182e6d933d577cda93` is the post-run `make_archive.py` hash; it is not the original `cross_face.py` run hash. This is a provenance distinction, not a mismatch in original-run identity.
- Archive bundle hash is `7be8ea23b3f1bbfddbe430ef97b63a135248072ae0bd9f60f41c1d79a295b536`; manifest SHA-256 `62e4f8948e6734a98f9660c211ccd7a86b22c4065243e1baf673d35612e4429c`; `archive_verify_receipt.json` records PASS for 15 cases. The portable verifier checks stored extrema, parameter membership, basic numeric arithmetic, statuses, and provenance. Its stated scope does not recalculate integrals or establish extrema over omitted raw rows; original raw rows are retained separately.

## MATLAB finite toy receipt and conclusion

The coordinator-provided `r047_primary_matlab_aggregation.m` and JSON receipt give exact rational values for a finite pinned-face aggregation toy for `n=2,…,7`. The code enumerates a finite three-state product and confirms its expansion identity. The receipt explicitly limits scope to exact finite-law aggregation: it does not evaluate the continuous R047 quadrature, certify a dimension-seven operator, or constitute `n=7` original-operator work. The `.m` file contains a historical private absolute receipt path and is retained unchanged; the portable archive verifier does not execute it.

The report's counts, extrema, resource figures, stabilization description, and stated limits agree with the registrations, summaries, monitors, retained raw cases, and archive metadata. No factual correction is indicated. R047 is a finite-grid 80-digit numerical diagnostic with reproducible retained inputs and post-run archival checks, not a proof of a global bound, sharp constant, or new `C_n` estimate.
