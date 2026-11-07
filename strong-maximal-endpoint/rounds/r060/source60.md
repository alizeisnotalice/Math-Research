# r060 source note — simultaneous source reoptimization

Assume the finite-measure/normalized-kernel hypotheses of r030, with
`0<q<cap`, so constant inputs supply Slater feasibility. Let `Phi` be the capped Orlicz integrand and let `R_u` be a family of
probability kernels.  This note gives the exact reduction that is relevant
when *the source itself is reoptimized* for several symmetric reception
constraints.  It is not the fixed-source response calculation of r059.

## Exact symmetric-pair reduction and defect formula

For an involution `sigma` preserving `d xi`, suppose
`R_- = R_+ o sigma`, and impose the two equal targets
`<H,R_+> >= q`, `<H,R_->=q`, with `0<q<cap`.  Assuming the usual capped
Fenchel strong duality, convex symmetrization gives the exact identity

\[
 P_{\{+,-\}}=P(\bar R),\qquad \bar R=(R_++R_-)/2,
\]

where `P(R)=inf{int Phi(H): <H,R> >=q}`.  Thus the *simultaneously
reoptimized* two-receiver source is the KKT source for `bar R`, rather than
the source of either receiver held fixed.

Write `lambda=e^ell` and, because the two individual kernels have the same
dual integral under the involution,
\[
 F(\ell)=e^\ell q-\int\Phi^*(e^\ell R_+)d\xi,
\]
\[
 D(\ell)=\int\left\{{\Phi^*(e^\ell R_+)+\Phi^*(e^\ell R_-)\over2}
 -\Phi^*(e^\ell\bar R)\right\}d\xi\ge0.
\]
Then, exactly (with no frozen-source approximation),
\[
 \boxed{P_{\{+,-\}}=\sup_\ell\{F(\ell)+D(\ell)\}.}\tag{60.1}
\]
If `ell_1` maximizes `F`, then
\[
 \boxed{P_{\{+,-\}}-P(R_+)\ge D(\ell_1).}\tag{60.2}
\]
This is a quantitative, nonnegative lower bound on the extra real source
cost caused by enforcing both receptions.  The remaining amount is exactly
the gain from reoptimizing the common multiplier in `F+D`; it has no fixed
sign representation from `D(ell_1)` alone.

For a finite symmetry group, or a compact symmetry group with normalized
Haar averaging, let `mu` be its invariant probability measure on a
transitive receiver orbit and assume the group preserves the base measure.
The same formula holds:
\[
 D_\mu(\ell)=\int\{\int\Phi^*(e^\ell R_u)d\mu(u)
 -\Phi^*(e^\ell\int R_ud\mu(u))\}d\xi\ge0,
\quad P_\mu=\sup_\ell(F+D_\mu).
\tag{60.3}
\]
This is a finite-amplitude covariance/Fenchel defect, not merely a
multi-multiplier Hessian.
For an arbitrary receiver set without a symmetry which makes all constraints
equal after symmetrizing `H`, (60.3) is only the averaged-constraint problem;
it is **not** equivalent to simultaneous feasibility.  That is the explicit
missing step for using a free continuous receiver measure as a many-constraint
upper or lower endpoint argument.

## Local link, and why it does not settle the natural window

For a smooth symmetric pair `R_±=R_0±epsilon S+O(epsilon^2)`, at a regular
KKT level,
\[
D(\ell)=\frac{\epsilon^2e^{2\ell}}2\,\mathfrak Q(S,S)+o(\epsilon^2),
\]
where `Q` includes both the lower contact measure and the soft inverse-curvature
term.  This recovers the *positive* piece behind the r055 pair calculation,
but (60.1) is finite-amplitude and simultaneously reoptimizes the source.
At a natural displacement of order `sqrt(N)`, that Taylor expansion is not
uniformly justified: a new estimate for the full finite defect `D_mu`,
including contact-set migration and the capped soft region, is still needed.

Consequently the mechanism obstruction is precise: convexity makes the
two-receiver source budget more expensive by (60.2), while a prospective
multi-cone construction needs enough union-volume growth to dominate this
defect.  No existing frozen-source `o(1)` action correction supplies a bound
on `D_mu` at `|u|~sqrt(N)`.

## Contact atoms in an unsplit discrete implementation

An attempted unsplit native Wolfram tensor calculation was discarded:
its binary choice at the lower KKT contact could not attain the target
mean. Such a grid replaces a contact layer by atoms; scalar bisection
alone does not solve that finite convex problem. No numerical source-cost
or defect claim from that attempt is retained.

The obstruction already occurs for one atom of unit mass, `R=1`, and
`0<q<1`: the exact optimum is `H=q` at `lambda=1`, whereas choosing only
`H=0` or `H=1` at contact is either infeasible or strictly suboptimal.
A discrete implementation must allow fractional contact mass; the
branch-split continuous evaluator used for the accepted many-receiver
search avoids substituting such an atomic binary problem.
