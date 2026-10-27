# C161 true second-count RQMC diagnostic — SOURCE PREPARED, UNRUN

No prepare script, random generator, kernel, special function, floating probe or old scientific baseline was run by Sol. Source identity, mathematical derivation, code writing, AST parsing and byte/hash bookkeeping alone constitute this packet. Terra must first receive separate main authorization for prepare-only, then a frozen-packet independent review plus main scientific release for the one finite batch.

## Fixed quantity and identity

The accepted original130 law, coefficients, body order and all21 cone weights remain unchanged:

    mu2 = kappa sum_(c=0)^20 alpha_c E[e2(k_h(B_c E))],
    kappa = 46080/pi^3, h=1/128, E_j iid Exp(1),
    e2(x_1,...,x_32) = sum_(i<j) x_i x_j.

The same six coordinates E drive all32 forms within a sample. The joint law is not replaced by independent body values, a pairwise surrogate, an orbit representative, or a density from149/151. The21 different cones are integrated separately with their original exact positive alpha values; all32 bodies and all496 pairs are included in each integrand. Source130 box_core/PROOF/results/inputs/cone_table/registration/source_manifest are read-only byte-pinned. Preparation rebinds every32x6 node to its original ray/rate and the body list (1,+/-1,...,+/-1), verifies positive alpha=det/product(rates), all21 ids and exact46080*sum alpha=R0. It references the accepted old determinant/fan proof rather than rerunning it.

No old343456 calculation, grid, kernel endpoint, MC sample, density, or baseline is replayed. The prepared source table is only model input. The optional question whether a fixed146 target near7.2e-8 is plausible does not alter h, weights, samples, sample count, family, or any existing computation.

## Stable floating kernel

Source130's complete closed expression for a>=0 is

    k_h(a) = 1/2 [erfc(u-s) - exp(2a) erfc(u+s)],
    u=a/(2s), s=sqrt(h).

By the definition erfcx(z)=exp(z^2)erfc(z) and

    2a - (u+s)^2 = -(u-s)^2,

this becomes exactly

    k_h(a) = 1/2 exp(-(u-s)^2)
                       [erfcx(u-s) - erfcx(u+s)].

This identity avoids the overflow-prone exp(2a) in the original floating expression. At a=0 it equals erf(s). All source coefficients are nonnegative and u>=0, so the left erfcx argument is bounded below by-s. Both erfcx terms and their difference are still ordinary floating evaluations with rounding, possible cancellation and tail underflow. No directed interval arithmetic, relative-error certificate, or proof of last-bit accuracy is claimed. Floating underflow to zero is allowed; invalid, divide-by-zero and overflow raise exceptions. The original proven0<=k<=p_h is used only as a numerical sanity check; there is no clipping, interval intersection, resampling or forced replacement of negative/nonfinite values.

All h, B and alpha entries are converted directly from original rational strings, not copied from display decimals. kappa_float is the direct conversion of the exact midpoint of130's saved rational kappa enclosure; that choice is explicitly numeric, not an outward bound. The exact symbolic kappa and enclosure remain in every final record. The p midpoint is a sanity reference only, never a proved floating bound. Source hash, version, conversion rules and relevant hex floats are recorded.

For each sample, set pair=prefix=0 and perform32 updates in the original body order:

    pair += prefix * k_i;
    prefix += k_i.

Induction gives pair=e2(k_1,...,k_i) over exact numbers. The floating implementation uses positive products/additions and avoids the cancellation in ((sum k)^2-sum k^2)/2. It does not make its sum rigorously outward-rounded.

## New finite randomized-QMC sample plan

The entire batch is fixed before science:8 replicates x21 cones x2^14 points of dimension6, hence168 jobs,2752512 points and16515072 uniform coordinates. Each job constructs scipy.stats.qmc.Sobol(d=6,scramble=True,bits=53,seed=seed) and calls random_base2(m=14) exactly once. The full chunk has16384 points. Independent-scramble terminology refers to separately seeded scrambled constructions, not an asserted theorem that deterministic hash seeds produce mathematically independent random variables. Each of168 distinct32-bit seeds is the first eight hexadecimal digits of SHA256 of its exact registered(rep,cone,source-table-hash) label. Preparation checks there is no seed collision and freezes the complete labels/digests/list. There is no seed search, discard, retry or adaptive extension.

Transform the resulting U in[0,1) by E=-log1p(-U). For an ideal uniform variable this gives Exp(1), coordinatewise. The floating Sobol net is a discrete approximation to that integral, as usual; we do not claim its finite empirical distribution equals the continuous law. In this transform U=0 maps to finite E=0, while the engine's53-bit construction excludes U=1. Do not shift, clip, replace or reject coordinates. Bounds and finiteness are checked. Each U matrix is hashed as little-endian float64 C-row-major bytes, and its min/max/zero-coordinate count, seed, algorithm, source and NumPy/SciPy versions are saved before kernel work. This avoids the singularity of using-log(U) at a zero coordinate.

Within each replicate r the estimate is

    z_r = kappa_float sum_c float(alpha_c) mean_(16384 points)[e2].

The final diagnostic is mean(z_0,...,z_7). Report all eight z_r, their empirical SD with denominator7, min/max and SD/sqrt(8), explicitly as descriptive dispersion. A QMC net's points are not treated as independent MC observations. No Student t interval, confidence level, statistical tail proof, rigorous lower/upper interval, or family exclusion is produced. Even if the estimate lies above the fixed146 target, this is only a reason to consider a separate future strict lower certificate; it cannot exclude the true value or authorize stopping151/expanding samples.

## Finite accounting and resource protocol

The all168-complete batch must exactly record:

- 168 Sobol job starts;
- 2752512 new sample points and16515072 uniform coordinates;
- 16515072 scalar log1p evaluations;
- 88080384 scalar kernel and exponential evaluations;
- 176160768 scalar erfcx evaluations;
- 88080384 prefix updates and88080384 pair-prefix updates.

These are vector-entry counts, not claims about implementation-internal CPU instructions or SciPy calls. The2^14-point engine call is one Sobol call per job; erfcx has two vector invocations per body, hence10752 vector invocations in total. The final sum includes496 pairs per sample algebraically, but no496-loop is executed. Counters increment just before work; partial/failed receipts therefore carry attempted or started counts. The worker records separate planned, Sobol-constructor-attempted, point-request-attempted, request-returned and complete-validated-sample states. Entering evaluate_job is not evidence of sampling. A version/API failure before a point request records zero requested samples; an attempted request without a validated return records actual sample counts null/UNKNOWN, never zero. Only fully completed jobs and aggregate equality establish the exact full counts. Hard process death makes uncheckpointed counts UNKNOWN.

At most15 spawned single-thread workers, nice>=19, standard BLAS/OpenMP thread variables1, no time/CPU ceiling, no automatic retry. Each point-axis array contains at most16384 rows; no full-batch point storage, pair tensor or grid is formed. Same-job U and E are16384x6; only one kernel vector is processed at a time over32 bodies. Each job gets exclusive output names and immutable source/sample/result receipts. A durable issued record precedes submit; ambiguous submit gaps remain owned/UNKNOWN. First failed job sets a stop event: no new job is issued, already running jobs may finish, not-yet-started jobs stop without sampling. The controller sets stop before cleanup on a body exception, then explicitly calls pool.shutdown(wait=True). Only an observed return from that call may set its reap observation to true; an interrupted or unobserved shutdown remains null/UNKNOWN. Child CPU is then an observed reaped subset unless shutdown was fully observed. Preserve all partial evidence; do not form a partial-cone mu2 estimate.

Save per-job PID/nice, CPU/wall/RSS, seed, U hash and finite counters. Controller saves start/end load, physical-page context, pool child CPU and immutable completion records; Terra separately records the actual terminal session exit/wall/reap receipt. Actual NumPy/SciPy/Python versions are captured by metadata at prepare-only, then checked in the scientific run. Unsupported API/version behavior stops, without a retry disguised as compatibility testing.

## Two distinct gates and claim boundary

Source deployment and prepare-only have a whitelist separate from scientific execution. prepare.py imports only the standard library, performs rational input identity and finite seed/registration bookkeeping, and makes no qmc/kernel calls. Its prepare_release.json must carry main approval and the exact source_freeze. It refuses overwrite/restart. After Terra freezes case_packet.json/registration.json, an independent review and main release must match all science requirements and choose workers<=15. Only then may run_prepared.py run once. No release file is included in this source packet.

This package supplies no numerical estimate yet. A successful future batch remains a numerical diagnostic of the true joint-law formula. It is neither a new strict mu2 lower/upper nor a new Cn/C6 bound. The established C6 upper remains unchanged. Existing130/156 strict evidence stays separate from the proposed descriptive estimate.

V2 is a new packet; original v1 manifest b28d475db71ac9299dfe6a43323dce61f351e1d075f31a67a337c9bc598c8787 stays unchanged. Only lifecycle/sample accounting and corresponding documentation changed. The numerical core from its np.errstate block through its return is byte-identical; plan.json, source_pins.json, prepare.py and common.py are byte-identical. Same168-job seed formula and finite sample plan. No prepare/science or exception-injection test was run in this revision.
