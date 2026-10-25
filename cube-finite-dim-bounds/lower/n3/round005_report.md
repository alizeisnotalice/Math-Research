# N3 round005 — certified satellite refinement and final 1000-evaluation batch

2026-10-25. **New lower bound beyond the full accurate round004 target T3:**

C3 >= 17940007113851367/10^16 = 1.7940007113851367.

The final saved union has a tight accurate interval in `round005c_certificate_result.json` and a still tighter independent interval in audit `007c_results.json`; its lower endpoint is strictly above the independently audited 005b complete upper endpoint, and above round004 T3's complete upper by more than 0.00002810698112149426. The simple displayed rational is accurately checked below that new lower endpoint. The independent audit007c entire final exact interval, rather than its shortened display, is the next comparison threshold. No optimality or upper-bound conclusion is made.

## Fixed background, exact construction, and immutable checkpoints

The background is Lebesgue and every witness is a centered axis-parallel cube with one common radius. The final rational construction is

mu = sum_{k∈Z^3} delta_{P*k} + epsilon sum_{k∈Z^3} delta_{P*(k+theta)},

with exact terminating-decimal rationals

- P=(1.630389466062437,1.072464196067522,1.129961851558576);
- theta=(0.613275060055078,0.659523766706894,0.053046037414171);
- epsilon=0.000731338850576.

Theta is the normalized phase, so physical displacement is P*theta coordinatewise. Unit primary mass and small satellite mass are both retained; mass per index pair is exactly 1+epsilon. All three periods, all three phases and the satellite mass changed during the final optimization. The final parameters were rationalized at 15 decimal places and all certification uses those saved exact values.

The original 170 witness keys remain saved throughout all three checkpoints. At the final geometry 158 boxes are nonempty; the other 12 original keys are retained with empty receiver boxes and contribute zero volume. `round005c_certificate_input.json` is self-contained; SHA-256 is `6d77724d8c4c7eaf416c9686e9c2cb6b041dfb2e11e0e0fe769a688fa7a64bc7`.

Original round005 and round005b inputs/results were preserved after their independent audit. Their input hashes are respectively `5ebe44bd14b85cb8483ac648a690f9226a07da35081eda848cd0ca38785f8cd7` and `8e26a3ad5e8d42f52ad3eb02855ac6721aabb67a75c50c7564e9751c82b13869`. They remain necessary intermediate checkpoints and seed/target evidence, not replacements for the final construction.

## Actual finite execution and stopping

Before starting, the exact target was read from master `results/bounds_ledger.json` and `audit/006-satellite_n3_results.json`; old round003/004/005 Python processes were absent. Neither optimization nor exact comparison used the short decimal display of T3 as its threshold.

|Stage|Actual optimization evaluations|Stopping reason|Full floating joint-union value|
|---|---:|---|---:|
|005 first checkpoint|2|gain >1e-8 above full round004 T3 upper; stop and certify|1.7939726890816543|
|005b, after independent checkpoint audit PASS|13|gain >1e-6 above complete 005 upper; stop and certify|1.793975867655425|
|005c final, after independent 005b audit PASS|1000|registered evaluation cap|1.7940007113851317|

Total actual optimization count is **1015**. No execution error, serialization recovery or replay occurred in this round. Early checkpoints interrupted their first Nelder–Mead simplex/search; the planned remaining starts were not run. Stages b and c were separately registered and explicitly authorized after the preceding independent exact audit. This certified continuation was not a silent continuation through an unverified candidate. The final 005c stage is the last stage of this round; no further optimization is running or authorized by its registration.

All objectives are complete finite three-dimensional receiver-box unions divided by **1+epsilon**. Sampling is used only to discover additional common-radius orbit-block keys, never as a sampled-volume estimate or objective. Stage 005 used 32768-point initial and independent survivor captures (seeds 20261025501/502). Stage b used 32768 initial and 65536 survivor centers (seeds 20261025511/512). Stage c used 32768 initial and 65536 final centers (seeds 20261025521/522). All six captures were completed; all found zero new keys beyond the saved 170 at those points. This negative finite capture result does not prove complete capture. Entry radii were bounded by 3.2 with offsets -4..4.

The final stage retained **every newly best value and parameter vector**, 188 records in the compact `round005c_best_trace.json`. It allowed gains below 1e-4 to accumulate through the full 1000 evaluations, per its registration. Its seed was exactly the saved 005b rational input. The fixed final-stage target came from the independent auditor's `audit/007b_results.json` full upper, frozen with source hash in `round005c_target.json`. The 1000-evaluation cap was reached without the >1e-4 early-stop threshold, and no convergence/local optimality is asserted. Initial/final witness capture and optimization took about 6.05 s for stage c.

Each evaluation regenerated all witness masses, critical radii and receiver geometries from the current seven parameters. Thus the same indexed family was not evaluated at stale radii or stale mass normalization. Original keys were retained even if they became empty. Recapture was performed at all selected checkpoints and final best; no new captured key was discarded before union evaluation.

## Exact box and union proof

A 12-integer key stores the inclusive lower triple for the primary orbit, lower triple for satellite, upper triple for primary, then upper triple for satellite. An empty orbit has lower=(0,0,0), upper=(-1,-1,-1). Let n_b=product_j max(0,h_bj-l_bj+1), M=n_0+epsilon*n_1, r=M^(1/3)/2, phi_0=0, phi_1=theta. Its receiver box in the full cell product[0,P_j] is the coordinate intersection of

[max_{b:n_b>0} P_j(h_bj+phi_bj)-r,
 min_{b:n_b>0} P_j(l_bj+phi_bj)+r].

For every point in this box the single cube of radius r contains all selected atoms. Unselected atoms can only add captured mass. The full joint receiver union is used; there is no product of marginal probabilities or inherited reflection shortcut.

`round005c_certify.py` reads the saved exact input without optimization. It uses the archived Fraction recursive union/root code in `certify.py`, explicitly a shared implementation. Cube roots are bracketed at denominator 10^24 with both cubed integer inequalities checked. Nested inner/outer boxes yield an exact ratio enclosure divided by 1+epsilon. An outward decimal display is

[1.7940007113851367947746180,
 1.7940007113851367947746211].

The exact rational lower is above audit007b's full upper by more than **0.00002484372971397893**. The accurate gap above the former round004 T3 full upper exceeds **0.00002810698112149426**. These comparisons are Fraction inequalities against saved full endpoints, not comparisons of displayed floats. The shorter rational 1.7940007113851367 is strictly between the fixed audit007b upper and this new lower endpoint.

`round005c_grid_check.py` is the second exact arithmetic path: SymPy integer_nthroot with integer inequalities, inward/outward 10^-12 endpoint rounding, then x-slabs with two-dimensional difference grids and Python-integer area/volume sums. It does not use the archived recursive union or its root helper. It gives a wider exact enclosure approximately [1.7940007113796217,1.7940007113905645], whose lower endpoint is also strictly above the complete independent audit007b upper. Both checkers use the same saved atom-to-box definition and input hash. The final discovery engine does not enter their verification.

Intermediate stage 005 has tight ratio about 1.79397268908165532947306 and accurate gain above T3 upper >8.46776400289595e-8. Stage b has tight ratio about 1.79397586765542281583630 and accurate gain above the full 005 author upper >3.17857376748636e-6; independent 007/007b audits separately verified their full target comparisons. The final source manifest pins the exact independently audited target and distinguishes it from the authors' slightly wider enclosures.

## Strict threshold, explicit index collar and finite L1 realization

The 170 saved keys have index_min=(-1,-2,-1), index_max=(2,2,2), hence **b=2**. The verifier computes these and checks every active orbit-block endpoint. Translate receiver cells by P*q, q∈[-N,N]^3. Each selected index k translates to k+q and lies within [-N-2,N+2]^3. Retain both original orbit representatives for every integer label in this block. Total atomic mass is exactly

(1+epsilon)(2N+5)^3.

Mass is paired by integer orbit label, not by physical-cell membership. The translated full cells have disjoint interiors, so received volume is (2N+1)^3 times the within-cell witness-union volume. Actual endpoint checks use N=0,1,7 and are saved in the exact result.

Nine saved radius records are critical equalities `(2*rlo)^3=M`. Closed critical boxes are not directly asserted inside a strict superlevel. For each nonempty retained receiver box, contract every face inward by tau. Its captured atoms lie at infinity-distance at most r-tau. Use the common-radius cube of radius r-tau/2 and smooth every retained atom by the same-mass nonnegative bump with support infinity-radius less than tau/2. Containment survives and (2r-tau)^3<M, hence the resulting L1 average is strictly above 1. One sufficiently small tau works for finitely many boxes. First N→infinity, then tau→0 recovers the critical joint ratio and proves C3 is at least the certified inner-union ratio. This uses the fixed Lebesgue background, actual finite-mass L1 functions and all-scale definition.

An explicit finite realization is also checked in `round005c_transfer_check.py`: take **N=10^7, tau=10^-12**, bump support infinity-radius tau/4. All periods are <2 and each grid-inner box has width >2tau. Shrinking a box costs at most 24tau in volume (six faces, each area at most 4); union shrink loss is at most the sum for 158 boxes. Therefore a rigorous finite weak-ratio lower is

(grid inner ratio − 24*158*tau/(1+epsilon)) * [(2N+1)/(2N+5)]^3.

The script computes this exactly and obtains approximately **1.7939996311904527**, still strictly above the complete audit007b upper. Thus the improvement is not dependent solely on an unquantified asymptotic transfer. Earlier checkpoints were independently audited before continuation; duplicate intermediate grid/finite scripts and MATLAB templates are not needed in the integration whitelist.

## MATLAB and remaining status

Stages 005 and 005b were actually checked in the shared local MATLAB before the main agent introduced serialized access. Their independently formed full grid-union values were respectively 1.7939726890811916 and 1.7939758676549566. These are floating cross-checks, not exact comparisons. Final stage c obeyed the updated rule: the main agent actually ran `round005c_matlab.m` serially and returned boxes=158, ratio=1.7940007113852767. This worker did not call shared MATLAB after the new rule. The floating result differs from the tight exact value by about 1.4e-13, consistent with floating root/summation effects and far below the accurate improvement margin. All scripts use distinct n3-prefixed variables and no clear/path mutation.

- Numerical discovery: the seven-parameter search, 188 retained final-stage best records, and six finite recaptures. No exhaustive or optimal-family claim.
- Local sufficient conditions: explicit atom blocks, one-radius containment intervals, b=2 collar, strict contraction and smoothing.
- Complete certification: two worker exact final finite-union enclosures, full-endpoint comparisons and explicit finite-L1 realization. Independent audit007c also passed on the same input hash: all 170 old keys were preserved, 158 boxes were nonempty, the other 12 were proved empty by outer bounds, 3009 atom occurrences/2653 nonempty endpoint checks were verified, and b=2/critical-level checks passed. That independently implemented certificate is separately owned by the main/audit agents and its full interval is the new target.
- Open: complete maximal-function superlevel set, additional capture shapes/scales, best satellite parameters and the optimal C3. No new upper bound or family exclusion.

No stage is still running. All old evidence and total ledgers were left untouched; no commits or pushes were made. The compact integration whitelist is in `round005_source_manifest.json`; it preserves necessary checkpoints, final source/input/certificates, best trace, target snapshot and run receipts, while omitting full logs/caches. Exact verification is read-only on the saved inputs. Reused helper/source hashes are explicitly listed rather than described as independent implementations.
