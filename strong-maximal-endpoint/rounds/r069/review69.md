# Independent final review of r069 forbidden-band capture and physical bridge

## Verdict

**PASS after one small but necessary proof repair in the capture corollary.**
The scalar forbidden-band asymptotic, the cap-two/cap-three deletion
argument, and the qualified fixed-action physical bridge are logically
consistent within their stated scope.  They do not establish that global
near-efficiency supplies the local-cost hypothesis, and they do not cover
receiver-dependent actions.  Both limitations are clearly retained.

## Scalar calibration and forbidden KKT problem

The proof correctly keeps the fixed-reference coefficient separate from the
general action-optimized one:

- fixed \(c=\mu\): \(\kappa=\beta=(3-\mu)/6\);
- deterministic shifted reference \(c_N=\mu-\sqrt{\mu}\rho/2\):
  \(\kappa=\alpha=(9-\mu)/24\).

The non-circular ordering is sound provided the cited r065 comparison is
read exactly as an \(L^1(P_c)\) comparison to the reference corner
profile.  It first confines the KKT cutoff to \(q_\circ+O(L)\), then
permits the local density and \(R=e^sf\) estimates, and only afterwards
identifies
\[
a_*=q_\circ+\kappa L+o(L).
\]
No general \(Q_N\) optimization coefficient is used in this step.

For the notched KKT problem, the identity
\[
Eh_{\lambda,I}=1-F(a)-P(S\in I)+S_I(a)
\]
is exact when \(a<\ell\).  The two comparison cutoffs
\(a_\pm=a_*-2KL\pm\varepsilon L\) do lie left of the band precisely because
\(\varepsilon<K-\kappa\).  Density flatness yields
\[
F(a_\pm)+P(S\in I)-F(a_*)
=\pm\varepsilon Lf(a_*)(1+o(1)),
\]
of order \(L^2/N\), while the soft/cap surplus error is only \(O(L/N)\).
This correctly brackets
\[
a_I=a_*-2KL+o(L).
\]
The retained left plateau has the claimed positive logarithmic length, and
the soft-plus-cap cost is \(O(e^{-a_I}/L)\).  Hence
\[
\log(J_{I,A}/J_{*,A})=2KL+o(L)
\]
is supported.  The cap-tail repair is adequate for every fixed \(A>1\):
the cap is reached only beyond a fixed positive multiple of \(N\), where
the stated Chernoff estimate is exponentially small.

## Capture calculation

For \(m=E[h1_I]\), \(m<1/3\), the rescaling
\[
g={h1_{I^c}\over1-m}
\]
satisfies \(0\le g\le3\) and \(Eg\ge1\).  The scaling inequality
\[
\Phi_N(ax)\le a^N\Phi_N(x),\qquad a\ge1,
\]
is valid across the kink as well as above it.  Therefore
\[
J_{I,3}\le J(g)\le(1-m)^{-N}J(h),
\]
and the cap-three forbidden optimum gives
\[
-N\log(1-m)\ge(2K-\varepsilon-o(1))L. \tag{R1}
\]

The displayed conclusion
\[
m\ge(2K-\varepsilon-o(1))L/N
\]
does follow, but one intermediate sentence is needed.  It is not valid to
replace \(-\log(1-m)\) by \(m\) solely from \(m<1/3\).  A correct short
argument is:

- if \(m\ge N^{-1/2}\), the claimed lower bound is automatic;
- otherwise \(m=o(1)\), so
  \(-\log(1-m)=m+O(m^2)=m+o(L/N)\), and (R1) gives the displayed bound.

This is a presentation repair, not a change in conclusion.

The band-cost inference is correct:
\[
m\le(\sup_I R)J_I(h),\qquad
\sup_I R=e^{q_\circ}N^K{L\over N}(1+o(1)).
\]
Together with
\(J_{*,2}=e^{-q_\circ}N^{-\kappa+o(1)}\), it gives
\[
J_I(h)\ge J_{*,2}N^{\kappa-K+o(1)}.
\]

## Physical bridge

For
\[
h(s)=K_c^sF_{u+s}(X,\sigma),
\qquad
j_u^c=\int e^{-s}K_c^s\Phi_N(F_{u+s})(X,\sigma)\,ds,
\]
conditional Jensen gives \(J(h)\le j_u^c\).  Thus the local hypothesis
\(j_u^c\le J_{*,2}N^\varepsilon\), together with receiver success, is
exactly sufficient to invoke scalar capture.  It is explicitly an extra
local assumption, not an implication of the global modular.

The step from band capture to
\[
e^{-u}\le Q_{\rm ref}N^{K-\kappa+o(1)}A_u^I
\]
is algebraically correct.  For a one-dimensional radial section,
\[
\int_{E_X}e^{-u}\,du\le\sup_{u\in E_X}e^{-u}
\]
is valid without halo monotonicity: a nonempty set is contained in
\([\inf E_X,\infty)\), whose exponential measure equals that supremum.
The empty section is zero.  Combined with the fixed-\(c\) r068 transport
estimate, this gives the qualified halo statement.  The note correctly
states that r028 is stronger and that this is a bridge-consistency check.

The \(1-\eta_N\) variant is also correct when
\(\eta_N=o(L/N)\): scaling costs
\((1-\eta_N)^{-N}=N^{o(1)}\), and a larger fixed repair cap handles the
new height.  The rogue small-scale dilation packet correctly shows why
global near-\(Q_{\rm ref}\) efficiency cannot imply capture at every
successful point while a common anchor remains fixed.

Finally, the prior-truncation comparison has the right direction:
the truncated KKT profile is feasible for the full law, and the dual
difference is bounded by \(A\lambda_{\rm tilde}d\).  Its stated scope
appropriately excludes floating/quadrature/tolerance certification.

## Required wording change

In the capture proof, insert the two-case argument following (R1) before
asserting the sharp first-order lower bound for \(m\).  With that addition,
I found no material gap in the stated qualified results.

Final integration note: the mass-bound repair is present as the exact inversion m>=1-exp[-(2K-epsilon-o(1))log(N)/N], giving the required first-order lower bound without assuming that m is small.
