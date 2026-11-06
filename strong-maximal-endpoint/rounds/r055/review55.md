# r055 independent review

## Local symmetric two-receiver identity

Let the two normalized kernels be related by coordinate swap,
`R_-epsilon=R_+epsilon^swap`.  Convex symmetrization of `H` preserves or
lowers the budget and turns the two target-one constraints into the single
mixture constraint for
\(\bar R_\epsilon=(R_\epsilon^++R_\epsilon^-)/2\).  Conversely a symmetric
mixture-feasible profile satisfies both original constraints.  Thus the two
infima agree.

For the contrast action `(c+h epsilon,c-h epsilon)` and translated anchors
`(a epsilon,-a epsilon)`, the mixture kernel is even in `epsilon`.  At a
regular KKT face its multiplier and profile have zero first variation.  The
pair second variation has no response term, whereas the single moving-kernel
budget has the usual nonnegative response.  With the stated sign convention,
\[
 P_{\rm pair}''=P_{\rm single}''+\lambda B_{\rm joint},\qquad B_{\rm joint}\ge0.
\]
Since both first derivatives vanish, division by the common value converts
this directly to the log-budget identity.  The r054 curvature supplies
\((\log P_{\rm single})''=Np_\mu h^2+o(Nh^2)\), provided `h` is fixed or
otherwise kept in the uniform local regime.

For positive `a epsilon`, the exact cone volume is
`V=2-exp(-2 a epsilon)`, and Wolfram 15 verified
\[
 \log V=2a\epsilon-4a^2\epsilon^2+O(\epsilon^3).
\]
Consequently the displayed score expansion with the extra
`lambda B/(2P)` term has the correct sign and factor `1/2`.

Required qualifications: use zero-extended block kernels on a fixed domain
and `b>=4`, so their required two derivatives have no moving-boundary delta
term; and state the result on a regular fixed KKT face.  It does not cover a
face transition without a one-sided/generalized derivative argument.

## Finite K-cone obstruction

For each individual receiver `u_j`, apply r026 to the translated restriction
`G(t)=H(u_j+t)`.  For balanced two blocks, cap two, and a single orthant
receiver this gives
\[
 J(H)\ge {e^{-\sum u_j}\over M_b Q_N},
 \qquad M_b=O(N^{-1/2}).
\]
Taking the strongest anchor and writing `w_j=e^{-sum u_j}` gives
\[
 {V\over J(H)}\le M_bQ_N E,
 \qquad E={V\over\max_jw_j}\le K
\]
for a union of `K` anchored orthants.  This is uniform in the fixed action
chosen for each receiver because r026 is uniform over deterministic ratio
vectors.  The conditions inherited from r026 remain essential: fixed block
count, balanced sizes, bounded/subexponential cap, nonnegative translated
orthant coordinates, and a **single-cone** receiver application.  It does
not apply to a full maximal halo.

For a decreasing piecewise-linear Pareto chain
\(u(v)=(d(v)+v,d(v)-v)\), put `w(v)=exp(-2d(v))`.  Segment integration and
one integration by parts give the exact formula
\[
 V_\Gamma={w(v_0)+w(v_1)\over2}+\int_{v_0}^{v_1}w(v)\,dv.
\]
Since `d` is linear on each segment, `w` has its maximum at a segment
endpoint.  Hence
\[
 E_\Gamma\le1+(v_1-v_0).
\]
This is the bound for the **continuous** chain; it is not bounded by the
number of its endpoint anchors.  For a finite set of `K` anchors,
`E_K<=K`.  If those finite anchors lie on the displayed Pareto chain, their
finite cone is contained in its continuous cone, so additionally
`E_K<=1+(v_1-v_0)`.  Therefore a finite-chain construction that reaches the
`Q_N` scale requires **both** `K` and angular span of order `sqrt(N)`;
gaining `N^gamma` requires both at least order `N^(1/2+gamma)`.  A genuinely
continuous chain has only the angular-span requirement.

## Verdict

The proposed identities and obstruction are valid with the listed
regularity/single-cone qualifications.  They do not claim a general
strong-maximal result or a multi-receiver lower bound.

## Audit of the final `local_geometry55.md`

Accepted after the final corrections.  The operator now correctly reads
\(\mathcal D_{a_0,h}=h(\partial_{c_1}-\partial_{c_2})
-a_0(\partial_{t_1}-\partial_{t_2})\), and the base kernel correctly
contains `e^(2d)`.  The response is unambiguously the r050 score-form
`B_{a_0,h}=mathfrak B(ell,ell)`, so (6) has the factor `lambda B/P`.

The strict-target issue is now handled correctly: all exact pair/curvature
identities and the first-order two-cone improvement apply at regular
`1<q<2`, but the r054 `Np_mu h^2+o(Nh^2)` asymptotic is invoked **only at
target one**.  Thus no unproved uniform-in-`q` curvature assertion remains.

Finally, the finite estimate explicitly uses the polygonal interpolation of
the finite vertices; each segment's maximum is at an endpoint, so the same
vertex maximum controls both the finite and polygonal cone volumes.  This
validates `E_finite<=min(K,1+Delta v)` in its stated polygonal scope.

## Signed-dual total-variation obstruction

This additional obstruction is valid in the same full-quadrant cap-two,
equal-two-block model.  Put
\[
 q(a)=a-\Phi_{N,[0,2]}^*(a),\qquad
 P_0=\nu\{H=0\}=E_\nu(H-1)_+.
\]
The second equality follows exactly from `E_nu H=1` and the KKT values
`H=0,1,` or `H>1`.  On `H>1`,
\[
 {q(a)\over a}=1-H+{\Phi_N(H)\over a},
 \qquad 0\le{\Phi_N(H)\over a}\le {4\over N},
\]
where the last bound holds on both soft and cap branches.  Hence
\[
 E_\nu(q/a)_-\ge P_0-{4\over N}.\tag{S1}
\]

The exact one-block integration-by-parts (Ward) identity is
\[
 {P\over\lambda}=E_\nu[(H-1)(1-k_1)].\tag{S2}
\]
It includes the lower jump: `d Phi(H)=a dH` there because the jump is from
0 to 1 at `a=1`.  With
`E(k_1-1)^2<=1/(b-2)` and
`E(H-1)^2<=2P_0`, Cauchy--Schwarz yields the unambiguous bound
\[
 {\lambda\over P}\ge\sqrt{\frac{b-2}{2P_0}}.\tag{S3}
\]
Since `q dxi/P=(lambda/P)(q/a)dnu`, its negative mass (and therefore its
total variation, up to an inessential factor two) is at least the product
of (S1) and (S3).

It remains to lower-bound `P_0`.  Profile stability gives `P_0=o(1)`, so
the KKT level `q_0=-log(lambda)` lies below a lower quartile of
`W=log R_b(T_1)+log R_b(T_2)`.  Uniform central local limits give
\[
 {W-Nm(c)+\log b\over\sqrt N}\Longrightarrow N(0,v_\mu),
\]
because `log f_{b,c}(T)+.5 log b=O_p(1)` on the central band.  Thus its
interquartile gap is at least `c sqrt(N)` for all large `N`.  On the upper
quartile, `W-q_0>=c sqrt(N)`, and the lower profile slope bound gives
\[
 (H-1)_+\ge\min\{(W-q_0-\log N)_+/N,1\}\ge {c\over\sqrt N}.
\]
Consequently `P_0>=c/sqrt(N)`.  Inserting this in (S1)--(S3) gives
\[
 \left\|{q(\lambda R)\,d\xi\over P}\right\|_{\rm TV}
 \ge cN^{1/4}.
\]
This is a real obstruction to treating the signed dual density as a
uniformly bounded-variation probability surrogate.  It does not affect the
restricted r054 curvature theorem and makes no assertion about `C_N`.
