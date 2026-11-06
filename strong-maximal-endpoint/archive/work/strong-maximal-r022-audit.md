# r022 independent theory audit — 2026-11-06

Read-only audit of `rounds/r022/report.md`; no repository/Git modifications.

## Verdict

Accepted. The proof establishes the stated sharp log coefficient for the
independent-proportion model (Q_n), with no transfer claimed to feedback or
general (C_n).

## Checks of the full chain

* The coarse inverse-modulus budget gives `F(s-2L)=O(sqrt(L/n))`; Chebyshev
  then excludes `D>=2d`, and Berry--Esseen locates `y0=M-s` at
  `sigma sqrt(nL)(1+o(1))`. This is sufficient before any local expansion.
* Under the negative tilt, the triangular Fourier local limit bounds give
  `f(s)/F(s)=lambda(1+O(L^-1/2))` and a global lower-side density envelope.
  The exact Gamma-plus-nonnegative-B identity gives
  `R(s)<=s f(s)/n`; it correctly has no omitted `B>s` term. It yields
  `F(s)=O(rho)`, negative excess `O(1)`, and then the essential sharp budget
  `F(s)<=y0/n+O(L/n)`.
* The all-vector Fourier CDF Edgeworth expansion has a genuine absolute
  `O(1/n)` remainder: uniform fourth moments control the local triangular
  expansion and the product Exp factors control the outer `/|t|` integral.
  At the test root this is `o(L^2/n)`, so the root inversion in §6 is valid.
* The deduction `D=O(L)` is noncircular: first (12) plus the coarse variance
  stability yields `D<=C L+C sqrt(DL)`, then (14) may be used.
* The global one-coordinate minimization of
  `H_lambda=D+lambda(v-mu)/2` gives the stated
  `-mu(mu-1)lambda^2/8` correction. A local Wolfram calculation independently
  returned `c_min=mu-mu lambda/2+O(lambda^2)` and that coefficient. Summing
  it produces exactly the `+(mu-1)L/8` variance compensation.
* Combining it with the Edgeworth `a=(3-mu)/6` term gives
  `(9-mu)/24`, matching the r019 lower witness.

The §8 Gaussian-initialization limitation and Gamma(2) splice observation are
properly labelled as mechanism scope, not used in the theorem. No blocker
found.
