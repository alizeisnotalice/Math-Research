# r059 independent audit

## Verdict: shifted continuation accepted; mesh needs two textual qualifications

Within their stated conditional scope (the r054/r057/r058 regularity,
compact-tilt, Ward, and finite-score estimates), I find no blocking gap in
the final shifted-response tail lemma, non-circular bootstrap, or Riccati
argument.  `mesh59.md` has the correct formulas and cost obstruction, but
its explanatory scaling after (4) must be corrected, and (5) needs an
explicit `Q>=2` hypothesis (automatic in the advertised asymptotic family).

## Shifted response and continuation: accepted

* From `lambda<=4b` and `R_c(t)<=2^b t^(b-1)/Gamma(b)`, the response support
  indeed gives `t1*t2>=eta*b^2`.  This transfers every source-axis inverse
  power to the opposite coordinate before the target Gamma upper tail is
  applied.  The proof correctly retains the density in the axis coordinate;
  it never uses a false global shifted pushforward-density supremum.
* The finite list `1,k_i,k_i*tilde{s}_j,k_i*(tilde{k}_j-1)` covers the
  response derivatives used here.  Target residuals are treated separately
  under the ordinary shifted target law, where the stated Chernoff/Gamma
  estimate applies.  Thus source `k` and target score/residual weights are
  not conflated.
* The bootstrap is non-circular: the provisional `B` is fixed first, the
  constants are independent of it for large `N`, and integrating the derived
  estimates strictly improves the provisional window.
* The final Riccati step is correct.  From `G(d(a,v),v,a)=1`,
  `d_aa=-G_aa/G_d-2(G_ad/G_d)d_a-(G_dd/G_d)d_a^2`.  On the stated action
  interval, `|d_a|=O(L^(3/2))`; hence the last two terms are `o(N)`, while
  `-G_aa/G_d~N`.  This establishes local convexity and the claimed
  `O(L^3/N)` fixed-versus-optimized action difference.

Minor clarity only: the axis paragraph should say explicitly that inverse
powers of `k` use the already stated lower bound `k>=(b-1)/t`, while positive
score factors use `k<=Cb/t`.

## Mesh: required correction

The exact crossing formula is
`eta*=(M-1)R/(1-MR)`.  Therefore, when `MR=o(1)`,
`eta*~MR`, **not** `R`.  This also follows directly from the suffix branch:
for small `eta`, its deficit is of size `eta/M`, so an order-`R` gap requires
`eta~MR`.  After undoing the outer normalization (whose physical outer
length is `M`), the physical separation is `M*eta*~M^2R`, not `MR`.

The formula for the adaptive benchmark `A*(eta)=1/(1+eta)` is consistent
with the stated adaptive comparison: each outer fibre may use its own inner
maximizing prefix of value one, while the nonadaptive common-prefix problem
has the suffix value `w/(w+eta)`.  Thus this is not a defect in (2)--(4),
only in the prose interpreting their scale.

For (5), `(1-1/Q)/M >= 1/(2M)` requires `Q>=2`.  It holds automatically
for the displayed asymptotic choice `Q=1+M*ceil((log N)^2)`, but should be
stated wherever (5) is asserted.

Apart from these two qualifications, the prefix identities, modular-cost
formula, endpoint reduction for outer escape, and recovery-cell lower bound
are algebraically consistent.  The note remains correctly limited to this
explicit conditional mesh family and makes no general `C_N` claim.

## Code receipt

`arc59.wl` imports `compute58.wl` under its library guard, freezes the source
profile, and differentiates the shifted target density.  Its four stored rows
converge; baseline `G_aa` and `G_d` checks agree with their independent Ward
references below `1e-8`.  This is appropriate numerical corroboration, not
a proof certificate.

Integration: both required mesh qualifications were corrected in the final note before commit; the source-axis score/Jacobian clarification was added.
