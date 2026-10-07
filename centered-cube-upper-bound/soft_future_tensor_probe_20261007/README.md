# Frozen full-tensor soft future-cap probe

2026-10-07. This directory owns only this new soft computation. It reads the complete old finite tensor input, lambda, hard winner and seed streams; it does not rerun old hard energy profiles or modify their data. Read the last sections of `checkpoint_20261004_general_flux_continuation.md`, the original soft/FIRST definitions in `概率接口阶段证明.tex`, and `auxiliary_exit_column_bridge_20261007.md`. The general weak(1,1) target remains unproved; this experiment supplies necessary-layer numerical diagnostics, without a new fee or dimensional fit.

Four receivers, two at n128 and two at n512, are selected **before evaluating any soft response**. For each old pair0/pair1 A-stream, choose the first saved time and first sampled source row passing hardband, the sufficient microbox count certificate, harddelay<=vstar and beta>=vstar. Full coordinates, old seeds, source component labels, hard winner and every coordinate arrival are in `frozen_points.json`. This deliberately conditional four-point catalogue does not estimate a source/cone population frequency or which sigma is generally common.

The source is exactly the old probability measure

    mu=sum_j p_j [ (64q_j)^(-1) sum_{k=0}^{64q_j-1} delta_{k/q_j} ]^tensor_n,
    q=(1,2,4), p=(1/4,1/2,1/4), W=1.

Resolution is shared by the entire source atom. Coincident component atoms carry summed masses. All components, boundaries and source masses remain. The common threshold is frozen q_first=lambda=64^(-n). No threshold or source parameter is optimized.

## Original soft kernel and full-grid evaluation

Use the original h=1_{[-1/2,1/2]}, w(t)=integral_0^1 exp(-|t|/v)dv, phi=h*w and psi_sigma=(1-sigma)h+sigma phi. For physical L in [1,2], normalized complete response is

    P_sigma,L mu(x)/lambda
      =sum_j p_j product_i [ (1-sigma)c_ji/(q_j L)
                            +sigma Sphi_ji/(q_j L) ].

Here c_ji counts the complete component coordinate grid in the closed receiving interval, and Sphi_ji sums phi((x_i-k/q_j)/L) over the complete finite grid. This is an exact tensor-mixture reduction, not an empirical source approximation. The original local quarter-unit source representation computes hard faces without subtracting large global coordinates.

For t=|x|, write g_a(v)=v exp(-a/v). The scalar phi integrand is 2v-g_{1/2-t}-g_{1/2+t} inside h, and g_{t-1/2}-g_{t+1/2} outside. Finite arithmetic-grid exponential sums are aggregated by geometric series, retaining the complete finite support. There is no spectral tail truncation: v covers [0,1], with integrand limit0 at v0.

For a>0, g_a''=a^2 v^(-3) exp(-a/v)>=0 and integral_0^1 g_a''=(1+a)exp(-a). At a0, g_a=v and its second derivative integral is0; using1 as an upper remains safe. Sum these positive curvature bounds for every grid term. The composite trapezoid Peano bound is

    |integral f - trapezoid_N f| <= (1/(8N^2)) integral |f''|.

Every coordinate scalar error is divided by q_j L and propagated through positive products and mixture sums. Three rounds use spectral N128/256/512. `validation.json` independently compares60 explicit-grid positive GL512 sums with geometric aggregation; all differences fall within the analytic scalar trapezoid remainder plus the stated numerical guard.

**Rounding remains unverified.** NumPy exp/log/geometric sums and products use binary64, with an additional documented1e-10 coordinate guard. That guard is not an outward-rounded interval proof. The Peano remainder is analytic in exact arithmetic, but none of the stored function enclosures or continuous cap statuses is a formal floating certificate.

## Hard jumps and continuous future coverage

All coordinate arrival radii of every full component are physical-L cuts, plus endpoints1/2 and a nested logL grid17/33/65. Later rounds retain earlier cuts and add up to64 targeted geometric midpoints per receiver. Every original hard jump has a saved strict-left/closed-right response trace. Final candidate-root jump sizes are approximately .1315,.1735,.4536,.2286; hard faces cannot be treated as a smooth L response.

The sigma grid combines k=n sigma in [0,32] with broad [0,1] nodes and larger dyadic k nodes. No monotonicity in sigma is assumed. Saved reverse cumulative maxima evaluate the entire **finite** future grid; a node above q excludes a proposed start numerically, while all finite nodes below q cannot prove continuous future cap.

A separate positive cell oracle covers full sigma/L rectangles. Each coordinate is affine in sigma, so the maximum of its two endpoint responses bounds the whole sigma interval; multiply these upper bounds and sum original component priors. Even decreasing h and phi give the physical-cell inequality

    sup_{l<=L<=u} P_sigma,L mu <= (u/l)^n P_sigma,u mu.

The closed upper endpoint retains all hard faces, even when several atoms enter at one L. Adaptive subdivision of sigma cells tries starts above the finite visible root by fixed dk/n gaps, dk=.5,1,2,4,8,16. A trial succeeds numerically only after exhausting all cells. Unresolved cells, early irreducible upper failure or the cell budget remain UNKNOWN. The sigma0=1 degenerate branch explicitly checks all L cells; `cap_endpoint_validation.json` tests that actual oracle branch with low/high synthetic coefficients. The bug fix and unaffected earlier trials are recorded in the round3 amendment; original outputs are retained.

The four final possible full last-crossing corridors are numerical only:

| n / pair | last clear finite-node violation lower | complete continuous numerical cap start upper |
|---|---:|---:|
| 128 / 0 | .6875 | .816771222 |
| 128 / 1 | .453125 | .515749621 |
| 512 / 0 | .15625 | .176822408 |
| 512 / 1 | .234375 | .248948026 |

These bracket the last-crossing candidate only conditional on the noncertified function enclosures. The displayed bisection roots .6917712/.4532496/.1611974/.2411355 are last-visible finite-L roots, **not interval boundaries or actual FIRST certification**. Bisection digits do not control the kernel integral. `root_error_brackets.json` additionally records fixed-L root uncertainty from the scalar integral bounds: final widths range1.15e-4 to2.78e-4. Actual soft physical winner Ls, possible hidden crossings, equality P=q and all GOOD/score/CP/GP/continuation gates remain uncertified.

## Weighted future average

At each displayed finite candidate pair, compute the original value

    K_sigma,L/q=alpha integral_0^1 v^(alpha-1)
      [P_{sigma+(1-sigma)v,L} mu/q] dv,
    alpha=ceil sqrt(n).

No new v quadrature is needed. The product is a degree-n Bernstein polynomial with positive normalized coefficients C_k. Its exact integration weights are

    weight_k=alpha/(n+alpha) product_{j=1}^{alpha-1} (k+j)/(n+j).

A positive coefficient recurrence and the scalar phi bounds produce the saved K lower/mid/upper values, subject to unverified binary64 rounding. Independent Gauss integration at a polynomial-exact order differs by at most4.94e-15 in `future_average_crosscheck.json`; root's separate saved review reports its own independent check.

Final finite-candidate K/q values are .432026/.273015/.00168826/.00303758. All scalar lower values exceed the theory's prechosen epsilon_n=log(n+2)^(-4); the n128 values also exceed1/16, and the n512 values do not. This is **not** a tested remaining low-future-average branch. K over the entire possible sigma/L candidate winner box is UNKNOWN. Actual FIRST and Ls are not certified, so high K at a finite candidate does not certify that actual original output is deleted. No theorem about the K column fee is proved here.

## Saved outputs and reproducibility

`summary.md` gives all three rounds and the numerical corridors. `completion_receipt.json` saves final scope and hashes. Each of12 NPZ files retains physical/sigma grids, scalar response lower/mid/upper, future-grid maxima and every hard-jump trace. Original coordinates/seeds, registrations/amendments, all per-point results and runtime receipts remain.

Round1/2/3 took71.44/59.98/112.98 seconds. All live run sessions terminated exit0. Earlier runner snapshots preserve the precise round1/round2 code hashes; the current runner includes the endpoint fix and runtime-only jump-distance caching. No finished round is restarted or overwritten.

Python used:

    /Users/zhengzhihao/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/bin/python3

Run each round separately only in a fresh copied output directory; existing round tags refuse overwrite. Scope is the original complete soft kernel and necessary shared future-cap layer. No actualgeom sample, weak-type bound, asymptotic exponent or new source-once spatial fee is claimed.
