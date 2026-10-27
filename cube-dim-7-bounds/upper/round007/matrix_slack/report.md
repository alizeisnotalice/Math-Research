# R007 matrix slack finite diagnostic

The sampled gap from full determinant to the eight-scalar-root bound is much larger than the subsequent scalar pair-AM–GM gap. This favors investigating matrix-level slack or four matrix pairs before investing in further scalar-moment refinement. The calculation is a finite pseudorandom diagnostic: none of the means or gaps is a deterministic expectation bound, a new strict I budget, a g certificate, or a C7 bound.

| observable | empirical mean | descriptive standard error |
|---|---:|---:|
| Z0: full matrix determinant | 0.0317397130 | 0.0000031640 |
| ZpairSum: sum of four pair determinants | 0.0538016220 | 0.0000033174 |
| Zpairmat: pair-determinant AM–GM | 0.0607357297 | 0.0000029314 |
| Z1: sum of eight scalar block roots | 0.0748028863 | 0.0000034011 |
| Z2: pairwise scalar AM–GM | 0.0810409271 | 0.0000029112 |
| Z3: eight scalar AM–GM | 0.0873047651 | 0.0000026117 |

Paired on each same shape, mean Z1−Z0≈0.04306317, versus Z2−Z1≈0.00623804, a ratio of about 6.9. Mean Z1−ZpairSum≈0.02100126 and Zpairmat−ZpairSum≈0.00693411; the sampled four-pair matrix budget is therefore a substantive diagnostic target. These are not justified subtractions from a certified upper bound. The incomparable difference Zpairmat−Z1 ranged from −0.03369 to +0.04069 in the sample; no total order is claimed. The Z2 and Z3 empirical means are consistent at ordinary descriptive sampling scales with the previously certified analytic expressions approximately 0.0810391403 and 0.0873027343, respectively; this is a numerical cross-check, not an error certificate.

The strict per-shape theory is Z0≤Z1≤Z2≤Z3, Z0≤ZpairSum≤Z1 and ZpairSum≤Zpairmat≤Z2, with each at most D/64. A fixed set of 15 independent PCG64 streams generated 2²⁰ original Gamma(.5,1) vectors each, normalized to Dirichlet(.5)^64, for 15,728,640 new shapes. Each stream has 1024 chunks and eight 2¹⁷-shape blocks. All six means, variances, extrema and all 15 paired signed differences, their stream/block summaries, and every task terminal are retained in results.json and summary.json. The displayed SEM assumes ideal iid sampling and ignores any separately uncertified floating arithmetic error; it does not bound rare shape regions deterministically.

The eight exact Hadamard block determinants were evaluated by their positive Cauchy–Binet formula, cross-checked against direct Cholesky on the first four already-counted shapes in each stream. The complete matrix and four matrix pairs used scaled Cholesky log determinants. A global Gamma max/min>2⁴⁰, small Cholesky diagonal ratio, nonfinite result or any raw hierarchy failure triggered same-input 128-decimal-digit recomputation, never clipping or resampling. There were 1997 such rows, below the shared hard cap 4096. Two initial float64 Z0≤ZpairSum violations exceeded the registered tolerance; both same-input high-precision checks passed. No final hierarchy violation was accepted. The 64 original Gamma values as exact floathex, their hashes, and stream/chunk/offset are retained for every fallback.

The batch used 24.220 s wall time and 153.416 s summed task CPU, about 6.33 core-equivalents over wall. Fifteen worker process snapshots totaled 608–690% CPU, with summed RSS 929728–1208672 KiB; summed RSS may include shared pages. The peak of any single worker was 95,764,480 bytes, and sampled system swap usage was 2332.38 MiB. The machine was shared with other work; these figures do not claim full utilization or a system-wide peak.

The frozen hierarchy proof and source registration are hash-linked. check.py validates the 15 task terminals and counts, all preserved fallback-input hashes, subblock and paired arithmetic, and output manifest. It does not replay the shapes, prove a floating-error allowance or certify an expectation. A rigorous improvement would require a separate deterministic lower certificate for matrix slack, or a certified matrix-pair moment/integration bound.
