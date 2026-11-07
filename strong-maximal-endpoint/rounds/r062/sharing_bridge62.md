# r062 — conditional sharing after the sharp response theorem

## Outcome

The fixed equal-block theorem now closes the former level-buffer gap. For
every coordinatewise nondecreasing cap A profile at critical modular cost
with one successful near-mu action, it gives uniform strict concavity on the
whole common window. Hence quadratic fibre regret and the weighted
barycentric sharing loss are derived under a fibrewise block
representation. The remaining bridge to uncentered rectangles is physical:
an endpoint-efficient label-preserving disintegration and a proof forcing
critical weighted variance. This does not bound C_N.

Put
\[
\rho_N=\sqrt{\log N/N},\qquad p_\mu=(\mu-1)/\mu,\qquad
I_N=[\mu-C\rho_N,\mu+C\rho_N].
\]

## 1. Uniform full-window curvature

For fixed r, equal blocks b=N/r, fixed cap A, and a coordinatewise
nondecreasing profile H, write \(\gamma_H(c)\) for its **common-mode**
response: all r block parameters equal the same scalar c. The sharpened
common-location result provides the scoped theorem:

> If
> \[
> J(H)\le\exp\{-N\mu+d_N+O(\log N)\}
> \]
> and one action \(c_*\in I_N\) succeeds, then the response is
> near-receiving throughout I_N and
> \[
> \gamma_H''(c)\le-\kappa_N,\qquad
> \kappa_N=p_\mu\log N(1+o(1)),\qquad c\in I_N. \tag{1}
> \]
> The estimate is uniform in H and does not use KKT.

The theorem follows from the sharp location response
\[
G_d(H,c)=(r+o(1)){\log N\over N},\quad
\|Z\|_2=O((\log N)^2),\quad
\|H-1\|_2=O(\rho_N^{1/2})
\]
and the common-mode density/Hessian reduction, with
\(m''(c)=-p_\mu+o(1)\). Its residual is o(log N). This is a common
parameter statement; no contrast direction is being used here.

The important new point is continuation. The rigidity proof extends to
\(EH\ge1-\eta\), \(\eta=o(\rho_N)\), and on that near-level component
\[
|\gamma_H'(c)|\le C\rho_N\log N. \tag{3}
\]
Starting at \(c_*\), variation over a window of width O(\rho_N) is
\[
O(\rho_N^2\log N)=O((\log N)^2/N)=o(N^{-1/2})=o(\rho_N). \tag{4}
\]
Thus all of I_N remains in the relaxed regime. No rescaling, modular budget
extension, or cap headroom is needed. This remains only a fixed
block/common-action theorem, not a theorem for physical rectangles.

## 2. Weighted barycentric sharing

Let \((Y,\nu)\) be the actual outer-fibre measure **normalized to a
probability measure**, including every Jacobian. (If a physical retained
piece has mass W, its unnormalized regret is W times the displayed result;
the lost mass must be charged in the endpoint accounting.) Suppose
\((y,c)\mapsto F_y(c)\) is jointly measurable, continuous in c, and has
the common-mode curvature (1) on \(I_N\). Define its constrained maximizer
\[
a_y=\underset{c\in I_N}{\operatorname{argmax}}\,F_y(c). \tag{5}
\]
Compactness and strict concavity give existence and uniqueness. Joint
measurability plus continuity makes \(a_y\) measurable (the measurable
maximum theorem, with a singleton argmax). If it is at an endpoint,
one-sided optimality still gives
\[
F_y'(a_y)(c-a_y)\le0\qquad(c\in I_N).
\]
Consequently strong concavity yields the pointwise, tie-free regret
\[
F_y(a_y)-F_y(c)\ge{\kappa_N\over2}(c-a_y)^2
\qquad(c\in I_N). \tag{6}
\]

With
\[
A=\int F_y(a_y)d\nu(y),\qquad
S=\sup_{c\in I_N}\int F_y(c)d\nu(y),
\]
integrating (6) and minimizing over the one common action proves
\[
\boxed{A-S\ge {\kappa_N\over2}\operatorname{Var}_\nu(a_y).} \tag{7}
\]
Varying fibre weights are already accounted for by the probability nu; no
arbitrary tie selection occurs.

If a physical representation supplies
\[
\operatorname{Var}_\nu(a_y)\ge v_0{\log N\over N}, \tag{8}
\]
then
\[
A-S\ge {p_\mu v_0\over2}{(\log N)^2\over N}(1+o(1)). \tag{9}
\]
For approximate maxima, an independently quantified fibre-regret error,
with weighted integral o(log^2(N)/N), is required. The main statement uses
exact unique maximizers.

## 3. Conditional vector sharing extension (awaiting audit)

The following is a conditional theorem based on `log_blocks62.md`; its
growing-block Hessian estimate still awaits independent audit. Let
\(L=\log N\), \(1\le r\le CL\), \(b=N/r\), and
\[
 \mathcal I_N=\{c\in\mathbb R^r:
 |c_i-\mu|\le K\sqrt{L/b}\ \hbox{for all }i\}.
\]
Assume that every jointly measurable fibre response
\((y,c)\mapsto F_y(c)\), arising from one successful capped critical
block profile, is continuous on \(\mathcal I_N\) and satisfies
\[
 \nabla_c^2F_y(c)\preceq
 -(p_\mu-o(1)){L\over r}I_r
 \qquad(c\in\mathcal I_N), \tag{10}
\]
uniformly in y. This is precisely the block conclusion proposed in
`log_blocks62.md`, not a proved geometric assertion here.

Let \(a_y\) be the unique constrained maximizer in \(\mathcal I_N\).
Measurability follows as in the scalar case. The box KKT inequality,
including faces and corners, is
\[
 \nabla F_y(a_y)\mathbin\cdot(c-a_y)\le0
 \qquad(c\in\mathcal I_N).
\]
Hence (10) gives
\[
 F_y(a_y)-F_y(c)\ge
 { (p_\mu-o(1))L\over2r}\|c-a_y\|_2^2.
\]
For the normalized actual outer-fibre probability measure nu, define
\(\operatorname{Cov}_\nu(a_y)\) in the usual vector sense. Integrating
and minimizing in the common vector c yields the conditional vector sharing
inequality
\[
 \boxed{\quad A-S\ge
 { (p_\mu-o(1))\log N\over2r}\,
 \operatorname{tr}\operatorname{Cov}_\nu(a_y).\quad} \tag{11}
\]
Thus the vector-critical condition
\[
 \operatorname{tr}\operatorname{Cov}_\nu(a_y)
 \ge v_0{r\log N\over N}
\]
again gives \((p_\mu v_0/2+o(1))(\log N)^2/N\). Until the analytic
estimate (10) receives its independent audit, this section is only a
conditional extension and must not be used as a proved block theorem.

## 4. The only remaining physical propositions

The general-rectangle bridge is now precisely:

> **P1 — efficient label-preserving fibre bridge.** Every
> endpoint-efficient actual input has a required weighted portion of its
> strict uncentered rectangle halo that disintegrates into fixed-r,
> balanced, near-mu block responses. The map preserves outer weights and
> every common action average \(F_y(c)\) on I_N, retains the critical
> modular budget and one success per fibre, and represents every competing
> genuine common rectangle in this action range. After normalizing the
> retained outer measure, the fibre response is jointly measurable in
> \((y,c)\) and continuous in c; the discarded mass is explicitly charged
> in the endpoint accounting.
>
> **P2 — variance forcing.** In the P1 representation, feedback-scale
> advantage forces
> \(\operatorname{Var}_\nu(a_y)\ge v_0\log N/N\) on the retained
> weighted mass, for fixed \(v_0>0\).

P1 includes common-action localization because it is a geometric part of
the label-preserving map, not a missing block analytic hypothesis. P1 and
P2 imply (9). They remain separate from final global endpoint accounting,
so no general claim for C_N follows here.

Label-preserving means retaining both physical side endpoints of every
rectangle while comparing all fibres with one common side/action. A
terminal-law or selected-action reduction cannot even form S in (7).

## 5. Audit of escape mechanisms

- r055's flat cube permits positive variance for an arbitrarily selected
  optimizer but zero common-action regret; one-sided outer intervals can
  also evade a mismatch found on a combined interval. It defeats
  unqualified variance selection, not the unique-maximizer theorem (7).

- r060 has a true critical log^2(N)/N common-action defect, low cost, and
  high conditional rank, yet bounded true strong-maximal score. Thus cap,
  cost, rank, and pointwise defect cannot imply P2; endpoint efficiency is
  essential.

- r030's iterated-maximal and origin-anchored reductions lose respectively
  common-side and translation labels, so neither can produce P1 directly.

These are mechanisms against unconditional bridges, not
endpoint-efficient counterexamples to P1/P2.

## 6. A concrete obstruction to naive conditional radialization

There is a further, elementary obstruction to one tempting proof of P1.
For one increment the density is
\[
 q_c(y)={e^{-y}(1+\mathbf 1_{\{y>c\}})\over1+e^{-c}}.
\]
Take two increments and condition their sum on \(S=s\). The conditional
density of \(Y_1\in[0,s]\), up to normalization, is
\[
 (1+\mathbf 1_{\{Y_1>c\}})
 (1+\mathbf 1_{\{s-Y_1>c\}}). \tag{12}
\]
Thus conditioning a fixed input on its radial sum generally creates an
action-dependent terminal profile.

For example, let \(s=3\) and use the fixed bounded input
\(h(y_1,y_2)=\mathbf1_{\{y_1<1/2\}}\). At \(c=1\), the denominator in
(12) integrates to \(8\), while the numerator integrates to \(1\), so
\[
 E_1[h\mid S=3]=\frac18.
\]
At \(c=2\), the corresponding denominator is \(5\), with the same
numerator \(1\), hence
\[
 E_2[h\mid S=3]=\frac15.
\]

The effect occurs inside the near-mu regime, rather than only at these
illustrative values. Set \(s=3\mu\), \(a=\mu/2\), and
\(h=\mathbf1_{\{Y_1<a\}}\). For c in a neighborhood of \(\mu\), one
has \(s>2c\) and \(a<c<s-a\). Formula (12) then has numerator \(2a\)
and denominator \(4s-4c\), giving
\[
 E_c[h\mid S=s]={a\over2(s-c)},\qquad
 {d\over dc}E_c[h\mid S=s]={a\over2(s-c)^2}\ne0. \tag{13}
\]
Therefore a reduction that radializes by conditioning on S cannot simply
declare the resulting terminal fixed while varying c. This refutes that
specific construction of P1, not the existence of some different
label-preserving, action-compatible P1 representation.

## Status

Within the capped fixed-block model, full-window quadratic fibre regret is
now proved. Only P1 plus P2 remain to connect it to general uncentered
rectangles; until then the result is conditional and does not upper-bound
C_N.
