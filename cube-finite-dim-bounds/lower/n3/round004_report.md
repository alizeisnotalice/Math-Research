# N3 round004 — a new certified small-satellite lower bound

2026-10-25. **The stated construction strictly improves the inherited accurate R3.** A simple exact lower bound is

C3 >= 17939726044040153/10^16 = 1.7939726044040153.

The saved 170-witness union has a tighter exact enclosure in `round004_certificate_result.json`. Its lower endpoint is strictly greater than both this simple rational and the original certificate's upper endpoint. The strict improvement over the original upper endpoint exceeds 0.00003937660924462117. No upper bound, optimality, or complete-superlevel-set claim is made.

## Exact construction and witness input

Write P=(1.630280147294,1.072431361926,1.129960333979), theta=(0.599580513532,0.664298207142,0.066908320568), epsilon=0.000561416422. Every terminating decimal here is an **exact rational**, not a rounded proof value. On the fixed Lebesgue background consider the periodic atomic input

mu = sum_{k∈Z^3} delta_{P*k} + epsilon sum_{k∈Z^3} delta_{P*(k+theta)}.

The common periods are P, not the old superperiod (2a1,a2,a3). Mass per integer orbit-index pair is 1+epsilon=500280708211/500000000000. Theta is a normalized phase; physical satellite displacement is P*theta coordinatewise. The positive satellite mass is genuinely small, with a unit-mass primary lattice. This differs from the old nearly equal two-orbit and three-layer families.

`round004_certificate_input.json` contains these exact rationals and all 170 witnesses; SHA-256 is `5f7c4e028a4b3d34d119410c24db9a569f96146427fdb60bd7eefed0e5515970`. Each key has 12 integers: primary lower triple, satellite lower triple, primary upper triple, satellite upper triple. A nonempty orbit selects its full inclusive integer block. Empty orbit marker is lower=(0,0,0), upper=(-1,-1,-1).

For each key, define n_b=product_j max(0,h_bj-l_bj+1), M=n_0+epsilon*n_1 and r=M^(1/3)/2. The receiver box consists of all x in the full cell product[0,P_j] satisfying

max_{b:n_b>0} P_j(h_bj+phi_bj)-r <= x_j <= min_{b:n_b>0} P_j(l_bj+phi_bj)+r,

with phi_0=0 and phi_1=theta. Every selected atom is inside the same centered cube of radius r. Extra unselected atoms only increase its captured mass. The full three-dimensional union of all 170 boxes is used; there is no marginal-product or reflection shortcut.

## Two exact union enclosures

The primary verifier `round004_certify.py` reads the saved exact input by default; it does not import or rerun the optimizer. It delegates Fraction root isolation and recursive joint-union arithmetic to the already archived `certify.py`, explicitly a shared implementation. Every mass is rational. Cube-root brackets use denominator 10^24 and both cubed integer inequalities are checked. Nested inner and outer boxes produce an exact ratio enclosure, divided by the full mass 1+epsilon. An outward decimal display is

[1.7939726044040153005134675,
 1.7939726044040153005134702].

The actual endpoints are the exact fractions in the certificate result and are the appropriate accurate next-round comparison target. The simple rational 1.7939726044040153 is strictly below the certified lower endpoint. All original comparisons read `baseline/original_certificates.json` directly, never rebuilding R3 from its displayed decimals.

A second checker, `round004_grid_check.py`, implements separate root and union paths. It uses SymPy integer_nthroot plus explicit integer inequalities for root isolation; rounds receiver endpoints inward/outward on the 10^-12 integer grid; and computes each integer union by x-slabs with two-dimensional difference grids and Python-integer area/volume sums. It does **not** call the old recursive union or its root helper. Its exact lower ratio is

1794979770081650568614032454922065759 /
1000561416422000000000000000000000000,

and its upper ratio is

1794979770086947012819258106627651413 /
1000561416422000000000000000000000000.

This wider enclosure is approximately [1.7939726044009219,1.7939726044062152]. Its lower endpoint is independently and strictly greater than the original accurate R3 upper endpoint. Both checkers verify the same saved input hash. Their common dependencies are exact rational arithmetic and the explicit atom-block-to-box definition; the discovery program does not enter either exact verification.

The inherited output fields `rigorous_gap_lower` and `strictly_below_original_R3` are retained only for compatibility. The former means original lower minus new upper, so its negative sign is expected; the latter is false. The authoritative improvement fields are `lower_above_original_R3_upper`, `rigorous_improvement_lower` and the checked simple rational lower.

## Strict superlevels, finite mass and general L1

Nine inner-radius bounds are exact critical equalities `(2*rlo)^3=M`. No closed critical receiver box is asserted to lie directly inside a strict superlevel set. For each retained nondegenerate receiver box, move every face inward by tau. Its selected atoms are then at infinity-distance at most r-tau. Use a witness cube of radius r-tau/2: it still contains all selected atoms, and its volume (2r-tau)^3 is strictly below M. Finitely many boxes permit a common sufficiently small positive tau; their union volume U_tau tends to the critical union volume U as tau decreases to zero.

The exact input keys give index_min=(-1,-2,-1), index_max=(2,2,2), hence **b=2**. The verifier computes these extrema and asserts every active selected index lies in [-b,b]^3. Receive physical cells P*(q+[0,1]^3), q∈[-N,N]^3. Translating one witness by P*q adds q to each orbit label. All required indices therefore belong to [-N-2,N+2]^3. Retain both original orbit representatives for every such integer label k. Total atomic input mass is exactly

(1+epsilon)*(2N+5)^3.

Group mass by orbit label k, not physical cell membership. In this construction both normalized phase representatives are already in the unit cell, but the mass proof requires only the integer labels. The actual check tests all active block endpoints at N=0,1,7: 912 checks passed. The received contracted union has volume (2N+1)^3 U_tau, up to shared cell boundaries of measure zero. Its weak ratio is at least

[U_tau/(1+epsilon)] * [(2N+1)/(2N+5)]^3.

Replace each retained atom by a same-mass nonnegative bump supported in infinity-radius strictly less than tau/2. The farthest support point remains inside the cube of radius r-tau/2; the strict volume/mass inequality remains valid. This is a finite-mass nonnegative L1 input on the fixed Lebesgue background. First let N tend to infinity, then tau decrease to zero. It follows that C3>=U/(1+epsilon), and thus C3 is at least each certified inner-union ratio. This establishes the new rational lower bound for general L1 and all positive scales, using only valid common-radius cubes.

## Executed search, stopping and calibration

The new interface started from the original algebraic spacing formula with u^3=2,v^3=3, not R3 decimals. Primary witnesses are the original selected count tuples translated across the full cell; satellite singleton translations are included. Mixed capture witnesses are regenerated with mass n_0+epsilon*n_1 and the new geometry whenever parameters change.

The registered scan used 32 scrambled Sobol phases (seed 20261025401) plus nine face/edge/center phases, at epsilon∈{1/256,1/64,1/16}: **123 cases**. Each used 4096 fixed discovery centers, seed 20261025402, to collect every passing orbit-block capture at entry radii <=2.8. The 660-key pool was then used in full finite-union optimization. Sampling discovers witness keys only; no sampled bad-set fraction or last-digit sampling gain is used as the objective or proof. Receiver-set ratios are always computed from the full finite joint box union and divided by 1+epsilon.

The planned cap was four 400-evaluation starts plus one 400-evaluation refinement. In fact only **800 optimization evaluations** ran: the first start ended at about 1.7938605603312479; the second ended at 1.7939726044037836 and triggered the registered immediate stop for exact certification. Both survivors were recaptured at 16384 independent centers, seed 20261025403. Both added zero keys to the pool. The remaining starts, final refinement and 32768-point capture were not run. Both completed starts stopped at their 400-evaluation cap, so neither is claimed locally optimal. Total discovery/optimization wall time was about 12.36 s. **There were no round004 execution errors, recovery replays or repeated optimization.** Earlier progress wording calling this the “first 400” survivor was corrected after the saved receipt confirmed it was the second start.

Optimization variables were all three log-period changes (bounded by .08), three free normalized satellite phases in [0,1], and log epsilon with epsilon∈[2^-12,2^-3]. The final actual discovery mass was about 0.000561416421958; certification uses the exact rational epsilon above. Bounds do not imply any exclusion outside this finite search. Discovery radii >2.8 were not searched; the resulting finite-witness proof does not require exhaustive capture.

The requested new-interface consistency calibration actually ran: with theta=0, both orbits select the same 37 primary blocks; multiplying all periods by (1+epsilon)^(1/3) exactly compensates for the combined coincident mass. At epsilon=1/256,1/64,1/16 the numerical joint ratios were 1.7939332277947715, 1.7939332277947724, 1.7939332277947702. Differences from the original coefficient-formula evaluation were <=2.3e-15. This checks the 1+epsilon denominator and counting of repeated-position mass; it is a floating interface regression, not a rerun of the old exact certificate.

MATLAB actually ran `round004_matlab.m`, independently built all 170 boxes and computed their joint union by a three-dimensional difference grid. It returned **1.7939726044036579**, within 3.6e-13 of the tight Fraction enclosure. That difference is consistent with floating root/summation error and is far below the exact improvement margin. MATLAB uses n3r4 variable prefixes and neither clears nor repaths the shared session. It is a floating cross-check; the two exact enclosures establish the strict comparison.

## Proof status, sources and integration

- Numerical discovery: 123 structured cases, two capped 400-evaluation starts, independent witness recapture and interface calibration. No exhaustive or optimal-family conclusion.
- Local sufficient conditions: the explicit 170 selected atom blocks, common-radius containment boxes, integer b=2, finite mass and contraction/smoothing construction.
- Complete certification: the two exact finite-union enclosures, strict comparison against the original R3 upper endpoint, and the simple rational 1.7939726044040153. This is a new C3 lower bound, subject to the stated analytic transfer above.
- Still open: the complete superlevel set of this periodic input, better satellite parameters, all other constructions and the optimal three-dimensional constant. No new upper bound or family exclusion is claimed.

No subsequent optimization is running. A later batch may use this saved exact construction as its seed and compare with the full exact new interval, but that work is not part of this batch.

Use the external Python venv. `round004_certify.py` and `round004_grid_check.py` verify the saved exact input without optimization. `round004_certify.py --build-input` is only an optional regeneration step after the search has regenerated its numerical survivor and convenience capture-key output. `round004_calibrate.py` regenerates the interface calibration. `round004_search.py` reproduces the finite search and immediate-excess stopping behavior. `round004_matlab.m` supplies the independent floating union check.

Only new `round004_` files listed in `round004_source_manifest.json` are integration artifacts, together with that manifest itself. The saved exact input includes all keys, so the duplicate round004_keys output is omitted from integration. Source hashes pin the shared Fraction verifier, numerical union helper and read-only master STATE/round003 report; they distinguish reused dependencies from the independent integer-grid checker. Temporary complete logs and caches are omitted. No old evidence, STATE, ledger or other project was edited. No commit or push was made.
