# r030: a geometry-consistent constraint between rectangles and scalar feedback

Date: 2026-11-06.  This note is intentionally about the *actual*
strong-rectangle geometry, not another fixed-ratio theorem.

## Exact rectangular-selector / zero-information lemma

Let \(f\geq0\) be bounded and compactly supported, and put
\[
 H_f=\{x:M_Sf(x)>1\}.
\]
For every \(x\in H_f\) choose an axis-parallel rectangle
\[
 R_x=\prod_{i=1}^n[a_i(x),b_i(x)]\ni x,\qquad
 {1\over |R_x|}\int_{R_x}f>1. \tag{G1}
\]
The choice can be made measurable after a null-set modification.  A
concrete construction is to enumerate rational-endpoint rectangles.  The
strict inequality in (G1), continuity of box integrals in their endpoints
for \(L^1_{\rm loc}\) inputs, and a small enlargement show that the
rational rectangles already cover \(H_f\); select the first one in the
enumeration.

Let \(\Theta=x\) be the receiver/rectangle label and let
\(U_1,\ldots,U_n\) be iid uniform variables on \((0,1)\), independent of
\(\Theta\).  Conditional on \(\Theta=x\), set
\[
 Y_i=a_i(x)+(b_i(x)-a_i(x))U_i.
\]
Then the rectangle witness is exactly
\[
 E[f(Y_1,\ldots,Y_n)\mid\Theta=x]
 ={1\over |R_x|}\int_{R_x}f>1. \tag{G2}
\]
In particular its side actions are all \(\sigma(\Theta)\)-measurable
*before any coordinate innovation is sampled*, and
\[
 {\cal L}(Y_1,\ldots,Y_n\mid\Theta)
 =\bigotimes_{i=1}^n{\rm Unif}[a_i(\Theta),b_i(\Theta)]. \tag{G3}
\]
Equivalently, with \({\cal F}_0=\sigma(\Theta)\), every physical
rectangle has the conditional zero-information property
\[
 I\big((a_i,b_i)_{i\le n};(U_i)_{i\le n}\mid{\cal F}_0\big)=0. \tag{G4}
\]

This is an exact, strictly narrower constraint than the arbitrary
predictable scalar policies in r011: a policy \(r_i=r_i(Z_1,\ldots,Z_{i-1})\)
generically has nonzero new information at stage \(i\), whereas a single
rectangle has already fixed every side at time zero.  The statement covers
all of \(M_S\), including nonanchored intervals and receiver-dependent
rectangles, because (G1)--(G3) are just its definition plus measurable
selection.  It is not restricted to r028's common anchor, reflected sides,
or fixed global ratio.

An equivalent analytic formulation is the shared-side constraint.  In two
dimensions a true rectangle contributes
\[
 \sup_{I_1\ni x_1,\ I_2\ni x_2}
 {1\over |I_2|}\int_{I_2}
 \left[{1\over |I_1|}\int_{I_1}f(y_1,y_2)\,dy_1\right]dy_2. \tag{G5}
\]
The same \(I_1\) must be used for every \(y_2\in I_2\).  The iterated
relaxation replaces this by
\[
 \sup_{I_2\ni x_2}{1\over |I_2|}\int_{I_2}
 \sup_{I_1\ni x_1}{1\over |I_1|}\int_{I_1}f(y_1,y_2)\,dy_1\,dy_2, \tag{G6}
\]
which permits the first side to depend on the later integration variable.
That is precisely the spurious information whose one-dimensional
rearrangement is represented by predictable feedback.

Thus a plausible intermediate upper model must retain an initial
receiver/rectangle label \(\Theta\) and the conditional-product
constraint (G3), or an equivalent multi-fibre shared-side constraint.  A
restriction only on the scalar state \(S_{i-1}\) cannot encode this:
the global rearrangement in r011 has discarded \(\Theta\) and the common
side information.

## Strictness is real, not an interpretation issue

For
\[
 f=10\,1_{[0,1]^2}+10\,1_{[1,2]^2},\qquad x=(3,3),
\]
the actual shared-side value is
\[
 M_Sf(3,3)=5/2,
\]
whereas the relaxed expression (G6), i.e. \(M_2M_1f(3,3)\), is
\[
 25/9>5/2. \tag{G7}
\]
The gap is \(5/18\).  It occurs because the first-coordinate maximizing
interval in (G6) changes between the \(y_2\in[0,1]\) and
\(y_2\in[1,2]\) fibres.  Therefore a sequential feedback / iterated-maximal
description cannot be treated as geometry-consistent merely because its
actions are predictable.

## A natural further simplification that is false

It is tempting to impose the zero-information constraint and then
normalize every physical interval to the reflected origin-anchored form
\([-rx,x]\), allowing \(r\) to vary freely with the receiver.  This still
does **not** cover general \(M_S\), already in one dimension.

Take
\[
 f=1_{\{1/2\le |y|\le1\}},\qquad x=5/4.
\]
The genuine interval \([1/2,5/4]\) gives
\[
 M_1f(5/4)\ge {1/2\over3/4}=2/3.
\]
In fact it is optimal.  For all \(0\le r\le1\),
\[
 {1\over(1+r)(5/4)}\int_{-5r/4}^{5/4}f
 \le {4\over9},
\]
with equality at \(r=4/5\).  Hence origin anchoring loses a strict
pointwise factor even when the ratio is allowed to depend on the receiver.
So an intermediate model covering \(C_N\) must retain a translation/sliding
side parameter in addition to any ratio parameter.

## Consequence for the r030 upper program

The useful concrete target is not a smaller class of ordinary scalar
Markov controls.  It is a **labelled open-loop product-kernel model**:
keep \(\Theta\), permit arbitrary translated side pairs
\((a_i(\Theta),b_i(\Theta))\), but require conditional product sampling
and shared sides as in (G3)/(G5).  It contains every actual strong-maximal
witness and excludes r011/r029's pathwise action adaptation.

What remains open is a global rearrangement/disintegration inequality that
preserves this label rather than collapsing it into the single decreasing
profile \(f^*\).  Without such an inequality one cannot turn (G3) into a
numerical bound between \(C_N\) and \(R_N\).  The two explicit examples
above show why both omitted structures—shared sides and translations—are
essential, so a proposed intermediate scalar model that removes either is
already invalid.

## Independent entropy obstruction: the scalar action kernel is not Fisher-smooth

This does not supply the desired geometry constraint, but it rules out a
natural proposed way of collapsing r011 feedback back to a fixed ratio.
For the one-step r011 increment law,
\[
q_c(z)={e^{-z}(1+1_{\{z\ge c\}})\over1+e^{-c}},\qquad z\ge0,
\]
the moving jump makes relative entropy linear, not quadratic, in an action
perturbation. For \(c>a\), direct integration gives
\[
D(\nu_c\Vert\nu_a)
=\log{1+e^{-a}\over1+e^{-c}}
-{\log2\,(e^{-a}-e^{-c})\over1+e^{-c}}, \tag{K1}
\]
whereas for \(c<a\),
\[
D(\nu_c\Vert\nu_a)
=\log{1+e^{-a}\over1+e^{-c}}
 +{2\log2\,(e^{-c}-e^{-a})\over1+e^{-c}}. \tag{K2}
\]
Writing \(p_a=e^{-a}/(1+e^{-a})\), the one-sided slopes at \(a\) are
\[
\partial_{c+}D(\nu_c\Vert\nu_a)|_{c=a}=p_a(1-\log2),\qquad
\partial_{c-}D(\nu_c\Vert\nu_a)|_{c=a}
=-p_a(2\log2-1). \tag{K3}
\]
Wolfram 15.0 actually differentiated (K1)--(K2); at \(a=1\) it returned
\(0.0825254334\) and \(-0.1038905545\), respectively. Thus there is a
strict cusp and no finite quadratic Fisher expansion at this parametrized
kernel family.

By entropy chain rule, comparison of a predictable policy to a fixed
action \(a\) pays
\[
 D(P_\pi\Vert P_a)=E_\pi\sum_i
 D(\nu_{c_i({\cal F}_{i-1})}\Vert\nu_a).
\]
An r029-type local action difference \(|c_i-a|\asymp i^{-1/2}\), on a
positive fraction of a time band, consequently costs
\(\sum i^{-1/2}\), not \(\sum i^{-1}\). It is far larger than logarithmic
at macroscopic horizons. Therefore a proof which asserts that feedback is
entropy-close to a fixed ratio by a quadratic Fisher-information bound is
invalid for these unsmoothed kernels.

This is only an obstruction to that scalar entropy-continuity method. It
neither proves that a feedback policy is physically realizable by
rectangles nor implies any separation involving \(C_N\).

## Final audit of r030_final report §§5--6

Date: 2026-11-06. **Accept.** For \(c>a\), put
\(r=e^{-a}\), \(s=e^{-c}\). Directly integrating the positive pieces of
\(q_c-q_a\) gives the exact total variation convention
\[
 {\rm TV}(\nu_c,\nu_a)
=(r-s)\left({2\over1+r}-{1\over1+s}\right)
=p_a(c-a)+O((c-a)^2). \tag{K4}
\]
This agrees with §5.  With the report's explicitly stated convention
\[
H^2(\nu_c,\nu_a)=\int(\sqrt{q_c}-\sqrt{q_a})^2\,dz
\]
(no factor \(1/2\)), the moving interval between the two jumps contributes
\[
p_a(\sqrt2-1)^2|c-a|+O((c-a)^2). \tag{K5}
\]
The smooth normalization changes contribute only quadratically, so there
is no omitted factor of two.

For deterministic \(c_i=a+\varepsilon/\sqrt i\), the right-hand KL
expansion gives
\[
\sum_{i=1}^N D(\nu_{c_i}\Vert\nu_a)
=p_a(1-\log2)\varepsilon\sum_{i=1}^Ni^{-1/2}
+O(\varepsilon^2\log N)
=2p_a(1-\log2)\varepsilon\sqrt N+O(\varepsilon^2\log N+1). \tag{K6}
\]
Thus the stated \(\sqrt N\) asymptotic is correct (with \(\varepsilon\)
fixed sufficiently small to remain in the local action neighborhood).
It correctly refutes only a quadratic-Fisher/entropy approximation
argument.

§6 accurately states the conditional-product/shared-side structure of
physical rectangles and explicitly says it has not yielded a scalar
intermediate upper bound or an independence constraint on r011 increments.
No implication about \(C_N\) is overstated.
