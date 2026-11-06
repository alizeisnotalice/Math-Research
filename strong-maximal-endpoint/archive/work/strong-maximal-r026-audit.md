# r026 balanced-block strict barrier audit — 2026-11-06

Audit based on the supplied theorem statement; formal r026 report was not yet
on disk. No repository/Git change.

## Verdict

The claimed inequality
\[
 J_m(H)^{-1}\le C n^{-(m-1)/2}Q_n
\]
is valid for fixed m, blocks `b_j>=kappa n`, bounded H with finite positive
budget, and `log A_n=o(n)`.

* The common tilt solving mean `s` exists uniquely. For
  `s in [n/2,3n]`, its rate `gamma` lies in `[1/3,2 mu]`; this gives uniform
  moments, a positive variance lower bound, and exponential characteristic
  damping even for unbounded c_i.
* Uniform triangular-array LLT gives total tilted density at its mean
  `>=c/sqrt(n)` and block density sup norms `<=C/sqrt(b_j)`. Including the
  determined final block in the product, the conditional `(m-1)` density is
  `O(n^{-(m-1)/2})`. Common exponential tilting leaves that conditional law
  unchanged.
* Jensen plus this conditional-density bound yields the stated factor in
  `J_m`. The original S leaves the central interval with exponential
  probability uniformly in the ratios. Scaling the central conditional
  profile by `(1-2 eps)^-1`, `eps=A_n e^{-kappa n}`, gives a strict scalar
  candidate; the Orlicz inflation is at most `a^n=1+o(1)`.
* The rearrangement step is legitimate because
  `e^s f_S(s)=E[e^B(s-B)_+^(n-1)]/Gamma(n)` is nondecreasing. Decreasing
  rearrangement in `x=e^{-s}` preserves scalar cost and can only increase the
  S-expectation, producing a nondecreasing Q-admissible terminal.

One drafting point: explicitly include the last block density in the displayed
conditional-density product, and state the usual measurability/finite-budget
assumptions. If the report also claims actual block-deep lower constructions,
state their separate bounded/truncation hypotheses; they are not needed for
the abstract barrier inequality above.

## Refined unbalanced statement

Also accepted. Each tilted block is `Gamma(b_j,rate gamma)` plus an independent
nonnegative Bernoulli shift, so its density is bounded by the Gamma density
supremum `C/sqrt(b_j)` (including `b_j=1`). Together with the total LLT lower
bound, this gives conditional density
`C_m sqrt(n)/sqrt(product b_j)` and hence
\[
 J_m^{-1}\le C_m\sqrt{n/\prod_jb_j}\,Q_n.
\]
This correctly reduces to the balanced factor and makes no false decay claim
when a block remains of fixed size.

## Formal report final audit

Accepted for theorem (1). The actual block-coordinate Jacobian is correct;
for a cone profile, exhaust the receiving cone by finite boxes before taking
their volume limit. §§4--5 correctly establish the common-tilt range, block
Gamma density bounds, total LLT lower bound, conditional-density factor,
central-tail loss, Orlicz scaling, and monotone weighted rearrangement.

Minor drafting clarification only: in §2, when discussing a root solution
with expectation exactly one under a hard height cap A, do not describe the
strictification as multiplying H if it can already attain A. Instead solve
the same capped convex problem with target `1+eta<A` (or increase its dual
parameter while retaining the cap) and send eta down to zero. The displayed
single-cone reduction itself assumes `E H(T)>1`, so theorem (1) is unaffected.
