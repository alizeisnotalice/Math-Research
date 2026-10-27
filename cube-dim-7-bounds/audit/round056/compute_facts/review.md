# R056 compute-facts review

**Facts disposition: ACCEPT.** This is a read-only audit of the actual finite numerical run, receipts, monitoring, and stated proof scope. No integration or MATLAB work was rerun; no author or main-clone files were changed. R057 is not included in this evidence.

## Frozen run identity and chunk coverage

Producer files are in `/Users/zhengzhihao/Desktop/T/cube-r030-compute/cube-dim-7-bounds/upper/round056/compute/`.

- `registration.json` SHA-256: `212ae10cde77a2efdcc48fbd33898a0f8b82fec60819126404a6c4d02c458da1`.
- Registered `r056_inverse.py` SHA-256: `0247dddec8ac652d44e3bfb35e8e3444585dd046d9c62f9712ac4273aa4db94f`; this matches the current source. The registered R051 proof source hash also matches the actual source file.
- `summary.json` SHA-256: `b0c5a21eb68703cc0726546cd46d9ba3448b6ad18c24e14ac38cb43650e2a6c2`; its registration hash matches `registration.json`.
- The summary lists 60 chunk files, all present, all status `completed`, each chunk SHA-256 matching the summary entry, with no missing failures. Their q-index ranges cover 0–127 exactly once for each of the 15 registered `d` values: `(0–31),(32–63),(64–95),(96–127)`.
- Per chunk, the stored receipt records 256 source states, 1,024 target combinations, 512 quadratures, 9,216 exact sampled column checks, PID, wall/CPU, peak RSS, maxima/errors and complete result rows. Summing those chunk fields supports 15,360 states, 61,440 target combinations, 30,720 quadratures, and 552,960 exact sampled column checks. Each chunk receipt is hashed in the summary; there is no separate chunk manifest.

The reported `g` correction is confirmed in both the frozen registration and actual code: `g=R−(κ−1)ρ=d−h+ζρ`, with `κρ=ρ+q−ζρ`. The registration hash matches the code used by the run, and the source computes `g=d−h+ζρ`; the brief contradictory initial task wording is explicitly corrected in the report. The producer separately confirmed that source/results were not modified after registration; the symbolic check is a post-run supplement.

## Resources

`summary.json` records 12.8032 s elapsed, 167.9585 s summed worker CPU, 15 distinct PIDs, and maximum per-worker RSS of 225,673,216 bytes (about 215.2 MiB). Summed worker CPU divided by wall time is about 13.12 core-equivalents as a batch average; this is distinct from instantaneous monitoring.

`monitor.jsonl` SHA-256 is `ecff6c34617b22d80c7a48fe72285dfb9b5adf64a728f8f96c41402a278d6ed8`; it contains three sampled snapshots. The highest observed aggregate CPU is 1446.8% (about 14.47 cores) at 10.64 s, with 2,969,248 KiB combined worker RSS (about 2.83 GiB). Earlier samples are 1035.3% at 0.54 s and 1440.0% at 5.58 s. The peak is a snapshot, not evidence of that utilization throughout the 12.80 s run or of memory exhaustion.

## Numerical scope and limits

- The run used 80-digit mpmath for per-source-state one-dimensional piecewise moment quadrature, comparing direct integrals to stated formulae. The stored variance is computed after normalizing the inverse-weighted law for each fixed source state. It is not a mixture variance over source states; mixing those states and using this variance as though it were a mixture variance would be unsupported.
- The 552,960 column checks use exact rational arithmetic at a finite set of points recorded per target. The elementary two-branch rationale is stated in the report, but finite point checks do not by themselves prove all geometric hypotheses or the continuous operator inequality.
- The largest recorded formula discrepancy and variance are finite-list observations, not certified global extrema or rigorous interval error bounds. The registered grid is finite and stops at the listed combinations. Passing every sampled check (and no observed value exceeding 4,096) does not prove a universal bound or a `C_n` estimate.
- `symbolic_identities.json` reports PASS for elementary identities and the sampled column argument; its own scope says geometric hypotheses require separate proof. It is post-run and not an independent author audit.
- This computation does not establish the R051/R056 flexible-budget theorem, a new `C_n` bound, a dimension-seven operator result, or lower-bound progress. The `n=ceil(64ρ/q)` values are parameters of this finite 1D diagnostic, not work on the original n-dimensional operator.

## Report/source provenance

`R056_REPORT.md` SHA-256 is `0a1e1f15804bc60cef915301e54b7541b5ea89475770618ee6b2a0f9311d1f6c`. The code and registration were timestamped at run setup; summary and monitor were created at run completion. `symbolic_identities.py/json` were created later and must be treated as supplemental evidence, not as part of the original registered batch. Report, registration, per-chunk receipts and summary agree on counts, corrected geometry, and scope. No factual correction is indicated.
