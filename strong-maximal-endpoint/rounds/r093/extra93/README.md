# r093 finite probe summary

This directory contains a reproducible streaming reducer and a compact summary. Raw JSONL stays in `/tmp`; it is not copied into the work tree.

## Model and probability audit

For the physical corner side ratio (r=e^{-c}), a sampled coordinate has the exact law

\[
Z_r\overset d=E+cB,\qquad E\sim\mathrm{Exp}(1),\quad B\sim\mathrm{Bernoulli}\!\left(\frac{r}{1+r}\right),\quad E\perp B.
\]

The probes use (p=1/(1+e^c)=r/(1+r)). They simulate (N-m) unit-weight coordinates in `base` and (m) coordinates in `part`, then compare `full=base+part` with `weighted=base+w*part`. This is an exact common-random-number coupling, with `full-weighted=(1-w)*part >= 0`. It is a subfamily of the fixed-side-vector model: all coordinates share the same side ratio and there are only two weight levels. It does not cover arbitrary heterogeneous side ratios or weights.

Each of the two empirical CDFs has (2{,}000{,}000) iid draws per job. For either law, DKW gives failure at most (2e^{-2S\eta^2}), with the source using (\eta=\sqrt{\log(4\cdot10^8)/(2S)}+2/S). The first term makes the two-law union at most (10^{-8}) per job; the added (2/S) covers the source's linear-quantile rank interpolation cushion. The DKW event is uniform over the CDF, so all four quantile probabilities are covered at once. The common-random-number dependence between the two laws does not affect the union bound. Across 23,324 weighted jobs, the family failure bound is at most (2.3324\cdot10^{-4}) (simultaneous coverage at least 99.976676%, under iid sampling and DKW). This only covers the finite input grid and its replicates.

The postprocessor normalizes both endpoints of each raw quantile-gap DKW interval by the same positive weight deficit \((1-w)m\) before comparing them to

\[
b_N=\log\frac{1+\tau_N}{1-\tau_N},\qquad \tau_N=\frac{\sqrt\mu\,z}{\sqrt N},\quad \bar\Phi(z)=\tau_N.
\]

## Critical-quantile results

Both weighted batches contain 23,324 jobs total, four replicates for every parameter model; the extension omitted the 680 jobs already computed in the first batch. The critical tail probability was separately retained from 0.02, 0.05, and 0.1 for each (N), and a unique critical row was found for every job.

Across all 23,324 critical rows:

- Minimum observed finite-difference `gap_per_weight_deficit`: **0.657667**, at (N=4,m=1,w=0.995,c=5.375,mathrm{rep}=2). Its corresponding (b_N=0.616086); the DKW interval after normalization is ([0,4.770738]).
- Smallest normalized DKW upper endpoint: **0.749504**, at (N=4,m=2,w=0.75,c=6,mathrm{rep}=0), with (b_N=0.616086). This is the closest interval upper endpoint to falling below (b_N), and still exceeds it by **0.133417**.
- Count with normalized DKW upper endpoint below (b_N): **0 of 23,324**.

| (N) | (b_N) | Minimum observed finite-difference ratio | Minimum normalized DKW upper |
|---:|---:|---:|---:|
| 4 | 0.616086 | 0.657667 | 0.749504 |
| 8 | 0.529514 | 0.723782 | 0.804287 |
| 16 | 0.445462 | 0.766418 | 0.844235 |
| 32 | 0.367447 | 0.777610 | 0.877869 |
| 64 | 0.297814 | 0.745108 | 0.905300 |
| 128 | 0.237689 | 0.758995 | 0.928349 |

The ten rows nearest the comparison boundary (sorted by normalized DKW upper minus (b_N)) are preserved in `summary93.json`, with ((N,m,w,c,mathrm{replicate})), seed, observed mean, interval, (b_N), margin, and SHA-256 of the original JSONL line.

These are finite Monte Carlo diagnostics for the selected two-level subfamily. They neither prove a uniform derivative lower bound nor establish a physical halo bound.

Across all four tail thresholds (0.02, 0.05, 0.1, and the dimension-specific critical probability), the minimum observed finite-difference sensitivity ratio is **0.219749** at (q=0.02,N=4,m=1,w=0.995,c=5.625,mathrm{rep}=0). Per-threshold minimum means and normalized interval endpoints are recorded in `summary93.json`. This ratio is not a derivative theorem; the broader 0.02/0.05/0.1 tails are auxiliary diagnostics and are not compared to the critical (b_N).

## Gaussian diagnostic and resources

The separate 768-row Gaussian batch covers (k=64,128,\ldots,4096) and (B\in\{1/16,1/8,1/4,1/2,1,2,4,8,16,32,64,128\}). Across 80/160 target precision pairs, the stored root residuals were below their targets (maximum residual/target ratio (5.59\times10^{-90})); compared values agreed at their saved precision, and all rank shifts were positive. `rank_shift_over_B` ranged from **1.006015** to **1.230298**. This is only a Gaussian root/rank-window diagnostic; it does not imply a physical onset dimension.

All three batches reported 12 actual workers, exit code 0, and no stop reason. Guard reports sampled macOS availability/RSS estimates, not hard memory reservations. The largest reported owned RSS was 1,883.9 MiB; minimum system-available memory was 9,994.2 MiB; each batch recorded a 1 GiB reserve. Summed elapsed batch wall time was about 498.8 seconds. The source did not record per-process CPU time. Detailed per-batch resource reports and SHA-256 hashes for scripts, model files, and all input JSONL/summary/resource files are in `summary93.json`.

## Reproducing the original local batches

Actual source copies: [Gaussian](original/fullcpu93.py), [initial weighted probe](original/weight_initial93.py), [extended weighted probe](original/weight_extended93.py). Copy each to its original `/tmp/strong_fullcpu93.py`, `/tmp/strong_weight_probe93.py` or `/tmp/strong_weight_probe93b.py` path before executing with the mathematics Python environment and `code/math_guard.py`. These source copies preserve the actual local paths and RNG specification. The reducer reads the original `/tmp` outputs; on another host adjust paths and record new hashes. The summary stores the original script and raw-input hashes, while the raw JSONL is intentionally outside the repository.
