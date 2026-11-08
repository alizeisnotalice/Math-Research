# Forbidden radial bands: scalar cost and qualified capture

## Scope and calibrated references

This is a scalar calculation at one reference action, not a moving-action
or general-\(C_N\) theorem.  A notched profile is not radially
nondecreasing, so it is not itself a monotone scalar witness.

Put \(L=\log N\), let \(S_N\) have density \(f\), CDF \(F\), and write
\(R(s)=e^sf(s)\), which is increasing for the Gamma--Bernoulli law.  Let
\(q_\circ\) be the uncorrected critical location and \(\kappa\) the
reference logarithmic coefficient.

For fixed \(c=\mu\), \(q_\circ=N\mu-d_N\) and
\(\kappa=\beta=(3-\mu)/6\), as in r019 Section 5 (11) and the
\(D=0,V=\mu N\) specialization of r022 Section 6 (12).  For
\[
 c_N=\mu-\frac{\sqrt\mu}{2}\sqrt{L/N},
\]
write \(c=\mu+t\sqrt{L/N}\).  Then
\[
 M=N\mu-\frac{p t^2}{2}L+o(L),\qquad
 V=\mu N+(\mu-1)t\sqrt{NL}+O(L),
\]
and r022 (12), (14) gives
\[
 \kappa(t)=\beta-\frac{p t^2}{2}
 -\frac{(\mu-1)t}{2\sqrt\mu}.
\]
At \(t=-\sqrt\mu/2\), this is
\(\alpha=(9-\mu)/24\).  The r019 shifted platform is the matching lower
reference, so this scalar reference matches \(Q_N\) to \(o(L)\), not
variable-action inputs.

## Non-circular scalar calibration

Let \(J_{*,A}\) be the unnotched scalar optimum with fixed cap \(A>1\),
and let \(a_*\) be its KKT cutoff \(\lambda R(a_*)=1\).  Capped r019
platform lower bounds and the specialized r022 upper bound first give
\[
 \log(1/J_{*,A})=q_\circ+\kappa L+o(L).                     \tag{1}
\]
This does not locate \(a_*\).  The one-block comparison in
r065/universal_action65, equation (9), gives
\[
 \|h_{*,A}-h_0\|_{P_c}=O(L^2/N).                            \tag{2}
\]
If \(a_*\notin q_\circ+[-CL,CL]\), large fixed \(C\) forces a unit profile
difference on a subwindow around \(q_\circ\) of probability
\(\asymp CL^2/N\), contradicting (2).  Hence \(a_*=q_\circ+O(L)\).
Only now local Cramér estimates imply
\[
 f(q_\circ+v)=f(q_\circ)(1+o(1)),\quad
 \partial_s\log R(q_\circ+v)=1+o(1)\ (|v|\le CL),\quad
 Nf(a_*)\asymp L.                                          \tag{3}
\]
The KKT estimate below yields \(J_{*,A}=e^{-a_*}(1+O_A(L^{-1}))\);
comparison with (1) proves \(a_*=q_\circ+\kappa L+o(L)\).

## Forbidden-band KKT calculation

Fix \(K>\kappa\) and
\[
 I=[\ell,r]=[q_\circ-KL,q_\circ+KL].
\]
Minimize \(J(h)=\int e^{-s}\Phi_N(h(s))\,ds\), with
\(\Phi_N(u)=u(1+\log^+u)^{N-1}\), subject to
\(0\le h\le A\), \(\mathbb Eh(S)=1\), and \(h=0\) on \(I\).
Outside \(I\), the pointwise KKT form is
\[
 h_{\lambda,I}=
 \begin{cases}
 0,&\lambda R<1,\\
 1,&1\le\lambda R\le N,\\
 (\Phi_N')^{-1}(\lambda R),&N<\lambda R<\Phi_N'(A),\\
 A,&\lambda R\ge\Phi_N'(A).
 \end{cases}
\]
The multiplier is new.  Let \(\lambda R(a)=1\).  For \(a<\ell\), set
\(S_I(a)=\mathbb E[(h_{\lambda,I}(S)-1)_+]\).  Then
\[
 \mathbb Eh_{\lambda,I}(S)=1-F(a)-\mathbb P(S\in I)+S_I(a). \tag{4}
\]
The unnotched identity is \(S_\varnothing(a_*)=F(a_*)\).

On the soft branch,
\[
 {\partial h\over\partial\log\lambda}
 ={h\over\partial_{\log h}\log\Phi_N'(h)}=O_A(N^{-1}).      \tag{5}
\]
For \(a=a_*+O_K(L)\), this gives
\[
 S_I(a)=F(a_*)+O_{A,K}(L/N).                                \tag{6}
\]
Indeed full soft surplus changes by that amount; on \(I\),
\(\lambda R=N^{O_K(1)}\), so the cap is absent and the deleted soft
surplus is smaller.

Choose \(0<\varepsilon<K-\kappa\) and
\(a_\pm=a_*-2KL\pm\varepsilon L\).  Both lie left of \(\ell\), so (4)
applies without a presumed cutoff location.  By (3),
\[
 F(a_\pm)+\mathbb P(S\in I)-F(a_*)
 =\pm\varepsilon Lf(a_*)(1+o(1)).                           \tag{7}
\]
This is \(\asymp L^2/N\), dominating (6), so the mean equation has
opposite signs at \(a_-\) and \(a_+\).  KKT monotonicity gives
\[
 a_I=a_*-2KL+o(L)<\ell.                                    \tag{8}
\]

If \(\lambda R(b)=N\), then \(b-a=L+o(L)\).  The left height-one plateau
has length at least \(\min\{1,K-\kappa\}L+o(L)\), contributing
\[
 e^{-a_I}\{1-N^{-\min(1,K-\kappa)+o(1)}\}.                  \tag{9}
\]
On the soft branch
\[
 \Phi_N(h)=\lambda Rh{1+\log h\over N+\log h}
 \le C_A\lambda Rh/N.
\]
The same bound holds at the cap, hence
\[
 J_{\rm soft+cap}\le C_A\lambda/N=O_{A,K}(e^{-a_I}/L).
\]
Therefore
\[
 J_{I,A}=e^{-a_I}(1+O_{A,K}(L^{-1})),\qquad
 \boxed{\log(J_{I,A}/J_{*,A})=2KL+o(L).}                    \tag{10}
\]

## Qualified scalar capture

Every fixed cap has the same \(\kappa\), by the repair below.  Hence
\[
 J_{I,3}=J_{*,2}N^{2K+o(1)}.                                \tag{11}
\]
If
\[
 0\le h\le2,\quad \mathbb Eh(S)\ge1,\quad
 J(h)\le J_{*,2}N^\varepsilon,\quad0\le\varepsilon<2K,
\]
and \(m=\mathbb E[h(S)\mathbf{1}_I(S)]\), then
\[
 \boxed{m\ge(2K-\varepsilon-o(1)){L\over N},\qquad
 J_I(h)\ge J_{*,2}N^{\kappa-K+o(1)}.}                       \tag{12}
\]
For \(m\ge1/3\) this is immediate.  Otherwise let
\[
 g={h\mathbf{1}_{I^c}\over1-m}.
\]
Then \(0\le g\le3\), \(\mathbb Eg\ge1\), and
\[
 \Phi_N(ax)\le a^N\Phi_N(x)\qquad(a\ge1,x\ge0),              \tag{13}
\]
by \(1+\log(ax)\le a\) below the kink and
\(1+\log(ax)\le a(1+\log x)\) above it.  Thus
\(J(g)\le(1-m)^{-N}J(h)\).  Applying (11) yields
\[
 N[-\log(1-m)]\ge(2K-\varepsilon-o(1))L.
\]
Exponentiating the preceding inequality gives directly
\[
 m\ge1-\exp\{-(2K-\varepsilon-o(1))L/N\}
 =(2K-\varepsilon-o(1))L/N.
\]
This exact inversion does not assume that m=o(1). Also \(\Phi_N(h)\ge h\), so
\[
 m\le(\sup_I R)J_I(h),\qquad
 \sup_I R=e^{q_\circ}N^K{L\over N}(1+o(1)).
\]
Since \(J_{*,2}=e^{-q_\circ}N^{-\kappa+o(1)}\), this proves (12).
Taking \(K=\kappa+\delta\) gives
\(J_I(h)\ge J_{*,2}N^{-\delta+o(1)}\).  This is fixed-reference scalar
capture only.

## Fixed-cap repair

For the r019 platform tail \(H_w(S-q)\), capping at any fixed \(A>1\)
changes it only on \(S\ge q+w+N\log A\).  Uniformly for \(c=\mu\) and
\(c_N=\mu-\frac{\sqrt\mu}{2}\sqrt{L/N}\), positive-moment Chernoff gives,
for \(q=q_\circ+O(L)\), \(w=O(L)\),
\[
 \mathbb E[H_w(S-q)\mathbf{1}_{\{S\ge q+w+N\log A\}}]
 \le e^{-\eta_A N}.                                        \tag{14}
\]
The exact r019 root derivative has size \(\gtrsim L/N\), so lowering \(q\)
by \(O(Ne^{-\eta_A N}/L)\) restores strict mean-one success and changes
\(\log J\) by \(o(1)\).  Thus every fixed cap has the same \(\kappa\),
justifying the cap-two/cap-three comparison in (11).
