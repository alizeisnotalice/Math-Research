# R002: finite-cache dangerous-point search and S3 pseudo-low diagnosis

Completed and stopped for primary review, 2026-10-27. This round establishes no strict global g floor, no global minimum, and no new C7 bound. The original generator, body64 weights, Lebesgue centered axis-parallel seven-cube, general L1 scope and all common positive radii are unchanged.

The principal result is a failure mode of the finite-shape search: all eight selected lowest optimizer points, both exact floathex coordinates and their rationalized versions, lie in S3 and have infinite true potential. Their low empirical values cannot refute the conditional target 1/179000. The eight fixed nearby points are exactly off S3, yet their finite-cache values barely change; these are numerical candidates showing that the current cache does not resolve the near-S3 tail, not credible certified low-G witnesses.

## Frozen scope and execution

Source commit: 0478c861dcbbfbbc8be472f7617ef4aa261f86eb. Independent clone: `/Users/zhengzhihao/Desktop/T/Math-Research-cube7-gprobe-r002`. Only the owned `cube-dim-7-bounds/upper/round002/g_probe/` path was written. The old clone and shape cache were read only. No commit, push, MATLAB, extra agent, cache reconstruction or new shape draw occurred.

Registration SHA256: `335a2f8528a6c150a7c5aacd0d409549c0ac48d40eddda826214d3af3850bad0`. Cache SHA256: `6dfa48bf8954a7d02270af5bda9c1d33d4fcaf0ba4f5e486589186b7f7c09552`; verified before and after the run. All evaluations use the same original 8192 shapes in eight original blocks of 1024, without dropping, clipping or resampling.

The 64 fixed unscrambled Sobol starts initialized L-BFGS-B on x=(1,y), y in [0,1]^6. The joint value and analytic gradient use the complete 8192-shape average. Optimization minimized log of this positive empirical value. Every start had a 256-call hard cap; none used more than 27 calls. SciPy reported success at all 64 starts: 53 projected-gradient stops and 11 relative-function-reduction stops. These are local finite-sample solver statuses, not mathematical optimality statements.

The first eight starts were completed and their S3 pseudo-low behavior was reported before the remaining registered 56 starts began. Persistent unique start handles, checkpoints and every optimization call are retained in `starts/000.json` through `starts/063.json`; no handle was re-dispatched.

| phase | starts / bundles | full 8192-shape point calls | additional high-precision subset verification calls | wall seconds | CPU seconds |
|---|---:|---:|---:|---:|---:|
| first stage | 8 starts | 130 | 0 | 0.402171542 | 0.257806 |
| second stage | 56 starts | 738 | 0 | 0.515214250 | 0.751727 |
| diagnostics | 8 bundles | 216 | 24 | 7.290010417 | 7.263902 |
| total | 64 starts, 8 bundles | 1084 | 24 | 8.207396209 | 8.273435 |

The 24 seven-shape high-precision checks were conservatively charged an additional point each: total budget charge 1108 of 32768. This leaves unused capacity, which is not authorization to expand the batch. Each diagnostic bundle contains original-point, rational-point and one fixed near-point evaluation, plus two full finite-difference stencils in six coordinates (27 full point calls per bundle).

The first stage used eight unique workers, peak eight concurrent handles; the second used fifteen unique workers, peak fifteen concurrent handles. BLAS/OpenMP limits were 1 in each worker. Measured worker-CPU / (configured workers * stage wall) was 8.01% and 9.73%, respectively. Jobs were short, so launch, cache load and coordination overhead limited occupancy; no high utilization is claimed. The final diagnostics ran in one process. Component timers include the recorded computations and exclude inter-stage human/agent communication gaps; their sum is not the entire chat elapsed time.

## Candidate results and exact singular-set classification

The table lists rationalized selected points, with free coordinates y only. Exact original optimizer floathex coordinates and every original rational representation are saved in the diagnostics and independent verification. Values are empirical averages of the original sample, not lower bounds.

| id | rational y | empirical G (1e-6 units) | empirical SEM (1e-6 units) | largest sample share of total | nearby off-S3 empirical G (1e-6 units) |
|---:|---|---:|---:|---:|---:|
| 0 | (1,1,1,0,1,1) | 1.515732098 | 0.087052557 | 2.239% | 1.515759366 |
| 1 | (1,1,0,1,1,387543/500000) | 1.764969960 | 0.111536057 | 3.869% | 1.765003781 |
| 2 | (1,1,0,1,773977/1000000,1) | 1.796060397 | 0.122420148 | 4.187% | 1.796094418 |
| 3 | (1,1,0,85821/125000,1,1) | 1.856024121 | 0.110447392 | 3.034% | 1.856056334 |
| 4 | (359309/500000,1,0,1,1,1) | 1.942381064 | 0.142174651 | 3.617% | 1.942409016 |
| 5 | (1,134111/250000,0,1,1,1) | 1.942467625 | 0.133089657 | 3.633% | 1.942499290 |
| 6 | (1,64069/1000000,1,1,1,44059/125000) | 1.964895202 | 0.122662250 | 3.246% | 1.964930867 |
| 7 | (1,0,1,543129/1000000,1,1) | 2.013913351 | 0.121397688 | 2.837% | 2.013945310 |

All eight rational points were classified by exact Fraction/Cramer arithmetic, with explicit body triples and coefficients retained. An independent checker directly verifies every witness and verifies that each exact floathex optimizer point lies in the same triple span. Thus there is no floating rank threshold behind the S3 classification.

For reference, body index b denotes (1, sign bits b_0,...,b_5), bit 1 meaning +1 and bit 0 meaning -1. The lowest rational point is already the average of body columns 55 and 63, hence lies in S2. The general S3 divergence statement follows by retaining a minimal independent r<=3 representation and conditioning its Cauchy coordinates within O(t) of their nonzero coefficients. Removing at most three body columns leaves a spanning remainder: its Gram matrix is 64Id minus at most three rank-one matrices of norm 7, so it remains at least 43Id. The remainder density has a local positive lower bound after t scaling. This gives p_t(x)>=const*t^(2r-7) and hence G(x)=infinity for r<=3. Negative representation coefficients present no issue for symmetric Cauchy densities.

The near-point perturbation is exactly y_k'=(1-7epsilon)y_k+epsilon*k, epsilon=1/1000003. Every near-point has seven distinct exact absolute coordinate values, and is therefore outside S3: three body directions can create at most four absolute coordinate values. Their full rationals are saved. Their finite-cache values are below the target as well, but the cache is continuous and smooth where all of its individual quadratic forms are positive; true G has singular tail behavior that this fixed finite sample can miss. Neither the singular-point pseudo-low nor its nearby continuation certifies a true dangerous point. The target's actual global validity remains unknown.

Eight original block means for the selected rational points, in units of 1e-6:

| id | block0 | block1 | block2 | block3 | block4 | block5 | block6 | block7 |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 0 | 1.454882 | 1.468209 | 1.384898 | 1.297212 | 1.514377 | 1.830073 | 1.704221 | 1.471985 |
| 1 | 1.495790 | 1.388666 | 1.507763 | 1.894710 | 1.509837 | 1.744128 | 2.520877 | 2.057989 |
| 2 | 1.424114 | 1.475494 | 1.353402 | 1.883867 | 1.561541 | 1.849835 | 2.723189 | 2.097041 |
| 3 | 1.606257 | 1.472251 | 1.508908 | 2.186305 | 1.829064 | 2.132726 | 2.476356 | 1.636326 |
| 4 | 1.435152 | 1.527829 | 1.445550 | 2.284575 | 1.835952 | 2.202791 | 2.757044 | 2.050156 |
| 5 | 1.577142 | 2.051273 | 1.531104 | 1.820454 | 1.922828 | 2.247268 | 2.109147 | 2.280526 |
| 6 | 1.611994 | 1.487191 | 1.569860 | 2.429759 | 2.236736 | 2.046175 | 1.694791 | 2.642656 |
| 7 | 1.668458 | 1.753636 | 1.752328 | 2.445709 | 2.273918 | 2.140681 | 2.064334 | 2.012244 |

Descriptive SEM and block stability do not bound the missing tail bias. In particular, these apparently moderate block fluctuations coexist with infinite true G at all eight central points.

## Arithmetic and gradient diagnostics

Every evaluated point retained all 8192 shapes. No nonpositive/nonfinite Q or value/gradient failure occurred. Original cache maximum condition estimate is approximately 6.24837e9. At each original, rational and near point, the seven rows with original condition>1e8 were recomputed at 80 decimal digits from the saved original P. The largest individual shape-contribution relative discrepancy was 1.16909154435e-6. These high-precision computations are diagnostics, not interval enclosures or a full expectation certificate.

Finite differences at h=1e-5 and h=3e-6 cross-check the analytic six-dimensional gradient, using central stencils or second-order one-sided stencils within the cube. Maximum absolute gradient discrepancy divided by G across selected candidates was 5.88447810909e-9 and 5.71435428506e-10, respectively. This supports the implementation of the finite-sample derivative. It says nothing about a true-potential lower floor or uniform Monte Carlo approximation near singular strata.

## Deliverables and reproduction

- `design.md`, `registration.json`: frozen problem, starts, caps, call accounting and pre-run hashes.
- `probe.py`: finite staged search, complete-sample joint value/gradient, rational rank checks and registered diagnostics. Existing stage receipts or start handles prevent silent repeats.
- `starts/*.json`: all 64 unique terminal handles, complete optimization traces, exact floathex points, solver messages, PIDs, timestamps and resource timers.
- `first8-receipt.json`, `rest56-receipt.json`, `diagnostics.json`, `receipt.json`: actual results, all block diagnostics, candidate coordinates, finite differences, high-precision checks and total accounting.
- `verify.py`, `verification.json`: independent saved-result accounting and exact witness checks; no kernel evaluations or new samples. Run with the registered local Python if needed.
- `reevaluate.py`: minimal read-only one-point reevaluation from the existing locked cache. It was not invoked during this experiment. For example, a future separately counted invocation may pass `1 1 1 1 0 1 1`; rationals or saved floathex strings are accepted. It prints one empirical value, gradient and block diagnostics, without modifying the cache.
- `sources-hash.json`: final manifest and read-only cache/proof provenance.

No phase restarts or extra experiments are pending. The finite batch is complete; additional tail-aware estimation, true g analysis or different candidate design requires a separately registered round.
