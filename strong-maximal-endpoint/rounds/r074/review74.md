# Independent audit: nested-rectangle action counterexample

## Verdict

**PASS.**  The nested-rectangle construction is a substantially simpler
counterexample to a *pointwise* actual-local-cost-to-fixed-reference
transfer.  It needs neither the posterior-tail lemma nor a radial baseline
comparison.  It is a single physical, full-sign-labelled source at the
fixed receiver \(x=(1,\ldots,1)\).  It does not refute an aggregate upper
bound, because it provides no common-source packing or global halo budget.

Let
\[
 R_c=[-e^{-c},1]^N,\qquad
 c_a=c_r+\Delta c,\qquad
 \Delta c={L^{1/4}\over\sqrt N},
\tag{1}
\]
where \(c_r\) is the optimized shifted reference.  For a fixed-cap scalar
minimizer \(h_a=h_{*,c_a,3/2}\), choose the standard wide central window
\(W_N\) on which the r065 radial-score bound is uniform and with
\(\mathbb P_{c_a}(S_{c_a}\notin W_N)\le2N^{-6}\).  Put
\[
 \eta_N=\mathbb E_{c_a}[h_a(S_{c_a}){\bf1}_{W_N^c}]
 \le3N^{-6},\qquad \vartheta_N={1+N^{-3}\over1-\eta_N},
\tag{2}
\]
and define in raw physical coordinates
\[
 f(y)=\vartheta_Nh_a\!\left(-\sum_{i=1}^N\log|y_i|\right)
 {\bf1}_{W_N}\!\left(-\sum_{i=1}^N\log|y_i|\right){\bf1}_{R_{c_a}}(y).
\tag{3}
\]
Give it value zero on coordinate hyperplanes.  This source is compact away
from coordinate axes, and \(\vartheta_N(3/2)<2\) for large \(N\).

## 1. Exact action law and local scalar cost

For a point chosen uniformly from \(R_c\), put
\(Z_i=-\log|Y_i|\).  Directly splitting the positive and negative pieces
of the interval gives
\[
 Z_i\ \stackrel{d}=\ E_i+cB_i,
 \qquad E_i\sim\operatorname{Exp}(1),\quad
 \mathbb P(B_i=1)={e^{-c}\over1+e^{-c}}={1\over1+e^c},
\tag{3}
\]
independently over \(i\).  Thus, with \(S_c=\sum_iZ_i\),
\[
 {1\over |R_c|}\int_{R_c}h(-\sum_i\log|y_i|)\,dy
 =\mathbb E_c h(S_c).                                    \tag{4}
\]
This is a direct Fubini identity in the original signed Cartesian
coordinates; no source label is changed when the action changes.

For \(c=c_a\), every point in the action's conditional angular fibre at
radial level \(s\) lies in \(R_{c_a}\), and the function in (3) is the
constant \(\vartheta_Nh_a(s){\bf1}_{W_N}(s)\) on that fibre.  Hence, in
the r069 local notation,
\[
 K_{c_a}^{s}f=\vartheta_Nh_a(s){\bf1}_{W_N}(s),\qquad
 K_{c_a}^{s}\Phi_N(f)=\Phi_N(\vartheta_Nh_a(s)){\bf1}_{W_N}(s).
\tag{5}
\]
Consequently
\[
 A_{R_{c_a}}f=\vartheta_N(1-\eta_N)=1+N^{-3},\qquad
 j_{c_a}\Phi_N(f)\le\vartheta_N^N J(h_a)=(1+o(1))J(h_a).
\tag{6}
\]
The latter is the *local cone* cost
\(\int e^{-s}K_{c_a}^s\Phi_N(f)\,ds\), not a claim about the full global
physical modular of (3).

Consequently this input has actual strict success and actual local cost
\[
 j_{c_a}\Phi_N(f)
 \le J_{*,2}(c_r)N^{o(1)}.                                \tag{7}
\]
Here the only scalar comparison needed is the already identified compact
action calculation: writing \(c=\mu+t\sqrt{L/N}\),
\(t_a=t_*+L^{-1/4}\) and \(\kappa'(t_*)=0\), the log-cost difference is
\(O(\sqrt L)+o(L)\).  The use of the \(3/2\)-cap minimizer on the left and
the cap-two reference on the right is essential; no same-action
\(1+o(1)\) cap assertion is made.

## 2. Exact reference collapse

Since \(c_a>c_r\), one has the literal set inclusion
\(R_{c_a}\subset R_{c_r}\).  The source is zero on the difference, so
there is no kernel comparison or posterior approximation to perform:
\[
 {1\over |R_{c_r}|}\int_{R_{c_r}}f(y)\,dy
 ={ |R_{c_a}|\over |R_{c_r}|}
   {1\over |R_{c_a}|}\int_{R_{c_a}}f(y)\,dy.
\tag{7}
\]
The second factor is the actual average from (4), including the harmless
strict and compact-tail repairs.  The volume ratio is exactly
\[
 { |R_{c_a}|\over |R_{c_r}|}
 =\left({1+e^{-c_a}\over1+e^{-c_r}}\right)^N,
\tag{8}
\]
and Taylor expansion of \(\log(1+e^{-c})\) yields
\[
 \log { |R_{c_a}|\over |R_{c_r}|}
 =-Np_{c_r}\Delta c+O(N\Delta c^2)
 =-p_\mu\sqrt N\,L^{1/4}+O(L^{3/4})\longrightarrow-\infty.\tag{9}
\]
Therefore the fixed-reference average in (7) tends to zero, not merely
below \(1-o(L/N)\).

## 3. The reference local cost also collapses

Let \(p_c(s)\) be the density of \(S_c\).  At a fixed \(s\), the
unnormalised signed angular slice of the *same physical source* (3) is
independent of whether it is subsequently tested in \(R_{c_a}\) or
\(R_{c_r}\).  Dividing the two slice identities by their rectangle volumes
gives the exact relation
\[
 K_{c_r}^s\Phi_N(f)
 =V_N\,{p_{c_a}(s)\over p_{c_r}(s)}K_{c_a}^s\Phi_N(f),
 \qquad V_N={|R_{c_a}|\over|R_{c_r}|}.                    \tag{10}
\]
On \(W_N\), r065 supplies the uniform radial-score estimate
\(\left|\partial_c\log p_c(s)\right|\le CL\).  Hence
\[
 {p_{c_a}(s)\over p_{c_r}(s)}
 =\exp\{O(\Delta c\,L)\}=1+o(1),\qquad s\in W_N,
\tag{11}
\]
and therefore
\[
 j_{c_r}\Phi_N(f)=V_N(1+o(1))j_{c_a}\Phi_N(f).            \tag{12}
\]
Thus even direct local-cost transfer across the two actions fails.  This is
still a one-receiver identity, not a global modular allocation.

## 4. Full labels and coordinate convention

Formula (3) already includes both signs in every coordinate.  A negative
coordinate has the raw domain \(-e^{-c}\le y_i<0\), which becomes
\(Z_i\ge c\); a positive coordinate becomes \(Z_i\ge0\).  Thus the
conditional Gamma--binomial mixture is only a disintegration of this one
physical rectangle.  In particular, it is legitimate that the source is
zero outside the **actual** rectangle \(R_{c_a}\): this is a fixed support
choice in raw Cartesian space, not a c-dependent redefinition when it is
tested by \(c_r\).  The reference test integrates exactly the same source
over the larger rectangle, which is precisely why (7) holds.

The window in (2) already supplies the compact repair: every nonzero
coordinate magnitude is bounded below by \(e^{-\sup W_N}\), on all sign
components.  This neither changes the inclusion nor affects (9).

## Scope

The example decisively rules out any theorem asserting, at each separately
successful receiver,
\[
 \text{strict actual success} +
 j_{c_a}\Phi_N(f)\le J_{*,2}(c_r)N^{o(1)}
 \quad\Longrightarrow\quad
 A_{R_{c_r}}f\ge1-o(L/N).
\]
It does not show that such rectangles can occupy a large aggregate halo
while sharing a global modular budget.  An aggregate upper proof can still
use action-indexed stopping, packing, or an anchor condition; those are not
addressed by this pointwise construction.

## Final integration precision

For the windowed profile, only the upper estimate j_ca<= (1+o(1))J_{*,3/2}(c_a) is needed. The scalar calibration is explicitly recorded with a compact-t uniform o(log N) remainder; its O(sqrt(log N)) quadratic contribution is not treated as a bound for that whole remainder.
