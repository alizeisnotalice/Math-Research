# r020 independent theory audit — 2026-11-06

Read-only audit of `rounds/r020/report.md`; no repository/Git modification.

## Verdict

Accepted. No Gaussian-comparison, barrier closure, or final probability-budget
blocker found.

## Critical checks

* In the one-step comparison, Taylor against `G~N(mu,mu)` yields the stated
  main term `D f + (v_r+D^2-mu)f'/2`. With
  `|v_r-mu|<=C_m sqrt(D)` and `R=|f'/f|`, the lower bound is
  `D/2-C_m sqrt(D)R/2>=-C_m^2R^2/8`; the square completion is exact.
  Weighted uniform moments and the explicit Gaussian density ratio control
  the unbounded Taylor remainder by `O(L f/tau)`.
* The right-shift density-ratio identity is correct. Its lower exponential
  bound is uniform after selecting `C0` after `K`, and gives shift gain at
  least `K L f/(4 tau)`. This dominates the comparison loss. Constants are
  noncircular: first `C_*`, then `K`, then `C0`, then `K0`; `H_N<=C_H L^2`
  has fixed `C_H`.
* Both outer bands close: left `W_k<=delta`; right uses the exponential tail
  of `Z_r` and a Gaussian split to make `1-EW_{k-1}<=delta`. The signed
  `V=W-epsilon-k delta` has `V_0<=1_{q>=0}` and preserves the recurrence,
  giving (7) for arbitrary predictable policies.
* Integrating (7) only to `T` creates `O(N^-3)` error; the positive-mgf tail
  handles the rest. The normal excess is `d+H_N` plus a lower-tail term, hence
  `d+O(L^2)`, without assuming adaptive total variance.
* The profile accounting is consistent. The base level is
  `P(S>s-2L)`, not `P(S<=s)`, and the middle inverse contributes the excess
  divided by `N`. The `b=sqrt(N)L^2` error is `d*b/N^2=O(L^(5/2)/N)`;
  low/high far-tail contributions are negligible. This proves (10).
* At the final test point, `K0>C_H+2` raises the threshold enough that the
  Gaussian left tail gains at least `L^2/(sigma sqrt(N))*phi(z_N)` which is
  of order `L^3/N`; this strictly dominates (10)'s `O(L^(5/2)/N)`.

Local Wolfram check reproduced the Gaussian density-ratio formula and the
square-completion inequality used in §3. No blocker found in the final
comparison or in the all-с initial/external-band treatment.
