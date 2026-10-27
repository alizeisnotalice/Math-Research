# C161 initial static review

Mathematical/sampling/claim-boundary review: PASS. Prepare-only can freeze the registered metadata without generating samples; final scientific gate awaits the concrete packet and the receipt qualifications below. Sol executed no RNG, NumPy/SciPy function, kernel, mathematical preflight or prepare script. Only text and SHA256/JSON metadata were inspected.

Source: `/Users/zhengzhihao/Desktop/T/Math-Research-cube6-sol-localmax-audit/cube-dim-6-bounds/upper/round161/true_beta2_rqmc_diagnostic`. Manifest `b28d475db71ac9299dfe6a43323dce61f351e1d075f31a67a337c9bc598c8787`; runner `ef72d1b40efd57a1ac73bf4e65bc17a0a99b0fed092344a5edffd674ca4125bf`.

## Fixed true-law diagnostic

The formula retains all21 original positive cone weights, all32 by6 body coefficients and the single original kappa factor. The same six E coordinates drive all32 body values within a row, so correlation of the physical linear forms is retained. Preparation verifies exact source body/ray/rate identities, alpha normalization, all IDs and source hash pins without regenerating the fan or calling an old functional. It inherits the already accepted determinant/fan proof. No149/151 density is consumed.

The stable identity is correct: for u=a/(2s), `2a-(u+s)^2=-(u-s)^2`. Multiplying the erfcx difference by `exp(-(u-s)^2)/2` is exactly the original whole signed expression. Its left argument is bounded below by -s, avoiding the exp(2a) overflow mechanism. Double-precision cancellation/underflow remain numerical errors; zeros from underflow are allowed, negative/nonfinite values fail, and no value is clipped or resampled. The saved p midpoint with tolerance is a sanity reference only.

The positive prefix recurrence gives every unordered pair product exactly once over exact numbers and avoids subtracting two close squared sums. Its floating result is not outward-rounded. Kappa is explicitly the float of the original accurate rational interval midpoint, alpha/body values are converted directly from rational strings, and no midpoint is promoted to a strict bound.

Exactly8 replicates by21 cones are registered, with one scrambled Sobol6 job of 2^14 points each, bits53. All168 unique deterministic seed labels/digests are frozen and collision-checked; there is no seed search. U=0 is legitimate and transforms to E=0. U<1 is checked before `-log1p(-U)`, with no clipping/shifting/rejection. This is a finite QMC approximation to the continuous iid Exp6 integral, not a claim that a discrete scrambled net has that exact empirical law or iid point observations. Distinct seeded scrambles are not asserted to be a proven mathematically independent random family.

The complete plan counts168 Sobol jobs,2752512 points,16515072 coordinates/log1p scalar entries,88080384 kernels/exponentials/prefix and pair-prefix updates,176160768 erfcx entries. These are vector-entry attempts/started counts on a partial path; only full job success and aggregate equality establish the complete exact counts. There is no496-loop or pair tensor. All168 results are mandatory before any estimate; a partial cone set is not substituted.

Replicate mean, denominator7 empirical SD, min/max and SD/sqrt8 are explicitly descriptive. No confidence level/t interval, rigorous uncertainty, true-value exclusion, beta2/Cn bound, or automatic expansion/stopping151 is claimed. Common floating bias is not measured by replicate dispersion, consistently with the noncertified interpretation.

## Resource gate qualifications

The single spawn pool has at most15 workers, nice19 and the five BLAS/OpenMP environment variables1 before scientific imports. Parent/worker source freezes, exclusive issued/job records and one-time markers prevent automatic retry. The issued record precedes submission, preserving ownership for an ambiguous submit gap. The default executor context waits for its original workers on ordinary completion, and worker failures set the shared stop event. Counts/proof progress after process loss must remain unknown.

Before final science release, amend the unconditional `all_original_pool_workers_reaped_by_context_exit=True` in the runner's final receipt. A KeyboardInterrupt can interrupt executor shutdown/join itself. Merely reaching `finally` does not prove complete reaping. At minimum record True only after a verified normal complete context exit and otherwise None/UNKNOWN; preserve original worker ownership. Setting stop inside the pool's exceptional path before shutdown/wait also prevents queued work from starting during a controller exception. No sampling/fixture is needed to review this control-flow fix.

The worker sets `started=True` before entering `evaluate_job`, whose imports/version checks and source conversions precede RNG construction. Thus its current `new_samples_started` field can be True even for a pre-RNG failure. Name that flag evaluation-entered/attempted, or derive sample-start evidence from the explicit Sobol counter/sample receipt and keep an interrupted construction unknown. Completed-job counters and the all168-success estimate are unaffected, but partial accounting must not label entry into a function as actual generated samples.

These are resource/evidence qualifications, not mathematical changes or permission for RNG/API testing. Source/version metadata and the fixed job list still need the separate Terra prepare-only receipt and subsequent exact-hash independent/main gate. No source, live generator, Git index or remote was modified by this reviewer.
