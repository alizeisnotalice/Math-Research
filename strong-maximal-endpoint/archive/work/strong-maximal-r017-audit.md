# r017 independent theory audit — 2026-11-06

Scope: read-only audit of `rounds/r017/report.md`; no repository/Git edit.

## Verdict

Accepted.  The passage from the fixed-ratio CDF result to the actual exponential
ramp score is mathematically justified and remains appropriately restricted to
that witness family.

## Critical checks

* At the soft root, `U=1+d/n+O(L/n)` is correct.  For
  `d=(1+eps)sigma sqrt(nL)`, the ordinary lower tail is `o(rho)`; for the
  `(1-eps)` comparison it dominates both `rho` and `n^-1/2`.  Since
  `0<=F-G<=E(q-S)_+/n<=O(n^-1/2)`, the two comparisons prove
  `d=sigma sqrt(nL)(1+o(1))`.
* The cited negative-mgf bound implies
  `E(q-S)_+ <= (nV/d) exp(-d^2/(2nV))=o(sqrt(nL))`.  Dividing by `n` gives
  `F-G=o(rho)`, so `A=1` indeed implies `F(q_n)~d/n~sigma rho`.
* The tilted density argument is noncircular.  First it obtains
  `f_N(a)/F_N(a)~lambda` from the exact one-sided tilted integral and the
  global Lipschitz bound.  Only then may the convolution Taylor formula yield
  `F_n/F_N=1+O(rho^2)`, and the analogous density formula gives
  `f_n/f_N=1+O(rho^2)`.  Thus `f(a)~L/n` follows legitimately.
* The derivatives
  `B'=-B/n-f`, `B''=B/n^2+f/n-f'`,
  `B'''=-B/n^3-f/n^2+f'/n-f''` are correct (independently symbolically
  differentiated).  With the r016 uniform density bounds, every term in the
  third-order Taylor remainder is bounded by
  `C[f(a)(lambda^2+n^-1)+n^-3]`; the exponential weights have uniform third
  moments.  This is `o(f(a)rho)`.
* Hence the numerator in (9) is
  `-f(a)[(c-mu)+mu lambda/2]+o(f(a)rho)`.  The otherwise material `1/n` in
  the denominator is smaller than `f(a)` by `1/L`, preserving the r016
  derivative coefficient.  Integrating uniformly over bounded `t` proves
  (2), and the (t=-sigma/2) score gain is `(mu-1)L/8+o(L)`.
* Formula (1) for the actual Orlicz denominator follows by `u=t-q` then
  `y=1+u/n`; it is c-independent.  The strict `q<q_n(c)` truncation/deep-layer
  argument therefore correctly realizes the boundary score only as a
  supremum, not as an invalid equality witness.

Actual local status: the existing Wolfram full-mixture 60/100-bit output is
now present for `n=256,1024`, with reported precision doubling stability; its
results agree with the displayed MATLAB values.  This is numerical support,
not used in the proof.

No blockers found in §2 coarse root, §3 tail density ratio, or §5 `B'''`
remainder control.
