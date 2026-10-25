# B008 finite registration: off-S2 integration

Owner: B Astra. Only this directory is writable. No stage/commit/push, no MATLAB, no delegated agents. Existing kernel baselines are read-only; no new kernel evaluation is planned.

Fixed task: body32 Cauchy frame in R6; q=6/5; I=int_[0,1]^5 G(1,y)^q dy. Current C6 upper bound remains 64. A local or scalar budget is not an improvement to C6.

Finite batch registered before numerical work:

- Four geometric distance layers delta=1/16, 1/8, 1/4, 1/2. These are distances in R6 to the union S2 of two-direction spans, restricted to the positive face.
- At most twelve scalar combinations: each of the four deltas with epsilon=1/1000, 1/100, 1/32. Compute explicit coarea/quotient bound constants and compare off-stratum short-time budget against the existing whole-face budget. Also evaluate the singular-neighborhood budget associated with each delta, from an analytically optimized one-parameter split (not a parameter search).
- Finite exact frame/combinatorial identities, dimension-dependent constants and derivative coefficients may be checked symbolically. These do not constitute spatial test nodes. No high-dimensional grid, Monte Carlo, or unbounded process.
- The explicit spatial-modulus constants and necessary uniform radius at the fixed accuracy target (1/11250)/10 are derived algebraically for the same four layers; they are not a new parameter search.

- No fresh kernel points. If an unexpected need for one arises, report to primary before changing this registration; six new points is the external ceiling, not a target.

Expected deliverable: self-contained proof of at-least-three-large-parameter event, iterated quotient/coarea bounds retaining both S1 and S2 distance, honest derivative limitations, a rigorous positive-face S2 tube-volume budget and all-face integration fallback. Quantify whether these constants can certify I^(5/6)<1/11250 (conditional on global g>=1/60000); otherwise state exact missing budget and minimum useful next step.

Inputs are frozen in input_manifest.json before the numerical script runs. Outputs include code, source/dependency hashes and versions, real runtime/PID/exit receipt. No whole-face certification is claimed unless coverage and every error term close.
