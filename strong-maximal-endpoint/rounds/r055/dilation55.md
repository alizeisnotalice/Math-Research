# r055 — dilation identity: signed-moment obstruction

## Scope

This note works on the full quadrant at a regular KKT face, and uses the r054
equal-two-block, single-receiver theorem
\((\log P)''_-=Np_\mu+o(N)\).  The contrast identity needed below is derived
here.  Nothing applies to a moving finite box, multiple receivers, a
nonsmooth KKT transition, or `C_N`.

## Exact sign of the dual density

For `1<h<2`, set `a(h)=Phi_N'(h)`.  Wolfram verified

\[
a(h)=(1+\log h)^{N-2}(N+\log h),\qquad
q(a(h))=(1+\log h)^{N-2}[N+\log h-(N-1)h]. \tag{1}
\]

The complete cap-two rule is
\[
q(a)=a\ (0\le a\le1),\quad q(a)=1\ (1\le a\le N),\quad
q(a)=\Phi_N(2)-a\ (a\ge\Phi_N'(2)), \tag{2}
\]
with the middle part parametrized by (1).  For `N>=3`, since
`d q(a(h))/dh=(1-h)a'(h)<0`, its unique zero after the plateau and before
the cap is given by \((N-1)h_*=N+\log h_*\): the bracket in (1) is one at
`h=1`, negative at `h=2`, and strictly decreasing.  Therefore

\[
q^+\le1,\qquad q\ge-a,\qquad q(a)\to-\infty. \tag{3}
\]

For `N=64`, Wolfram gives
\(h_*=1.01612695712943292018\ldots\) and
\(a_*=171.25897446122018941\ldots\).  Thus the signed region begins very
close to the soft threshold; the cap branch remains negative.

## Self-contained contrast identity

Let `N=2b`, pull back by `t_i=c_i x_i`, and set
\(R_{b,c}(cx)=c^{b-1}S_b(x)/(\Gamma(b)Z(c)^b)\).  Along
`c_1=c+theta,c_2=c-theta`, write `D_x=x_1-x_2`.  Direct logarithmic
differentiation of the pulled-back exponential and Gamma--Bernoulli factors
gives at zero
\[
{\xi'\over\xi}={\nu'\over\nu}=-D_x,\quad
{\xi''\over\xi}=D_x^2-{2\over c^2},\quad
{\nu''\over\nu}=D_x^2-{Nv(c)\over c^2}. \tag{4}
\]
The same pullback gives `(log R)'=0`.  Swap symmetry at a regular face gives
`lambda'=a'=H'=0`; hence the envelope derivative has no response term:
\[
P''=\int\Phi(H)\xi''-\lambda\int H\nu''. \tag{5}
\]
Fenchel equality gives `Phi(H)-aH=-Phi^*(a)`, while the mean-one constraint
and (4) give \(\int D_x^2a\,d\xi=\lambda Nv(c)/c^2\).  Substitution in
(5), using \(q=a-\Phi^*(a)\), yields
\[
 (\log P)''_-={1\over c^2}\left\{{\int(t_1-t_2)^2q\,d\xi\over P}-2\right\}. \tag{6}
\]
This derivation is full-quadrant only; a moving finite box has boundary flux.

## Consequence of positive Np curvature

Put \(D=t_1-t_2\) and \(M_N=\int D^2q\,d\xi\).  Equation (6) and the r054
curvature theorem give
\[
M_N=P\{2+c^2Np_\mu+o(N)\}. \tag{7}
\]
In the original `t` coordinates, `dξ=e^{-t_1-t_2}dt_1dt_2` is the product
of two unit exponential laws.  Wolfram checked \(\int D^2d\xi=2\).  Hence,
by (3),
\[
M_N\le\int D^2q^+d\xi\le2,\qquad
P\le {2\over c^2Np_\mu+o(N)}=O(N^{-1}). \tag{8}
\]
This is a real but one-sided normalization consequence.  For fixed witness
volume it yields a *lower* score restriction, not an upper score bound.

The exact tail capacity is
\[
\int_{|D|>L}D^2d\xi=(L^2+2L+2)e^{-L}. \tag{9}
\]
Thus `M_N>r(L)` implies only
\[
\xi\{q>0\}\ge {M_N-r(L)\over L^2},\qquad r(L)=(L^2+2L+2)e^{-L}. \tag{10}
\]
It is not a fixed positive probability because its scale still contains the
unknown `NP`.

## Total-variation obstruction

There is a stronger consequence under the profile hypotheses of
r054.  Fenchel equality gives exactly
\[
{q(a)\over a}=1-H+{\Phi_N(H)\over a}. \tag{11}
\]
On the soft/cap part `H>1`,
\(\Phi_N(H)/a\le4/N\).  If
\(P_0=E_\nu(H-1)_+\), then
\[
\int q^-d\xi\ge\lambda(P_0-4/N). \tag{12}
\]
Write H=0 below the KKT level and H>=1 above it. Since E_nu H=1, P_0=nu(H=0) and E|H-1|^2<=2P_0. The Ward identity P/lambda=E[(H-1)(1-k_1)], E k_1=1, and E(k_1-1)^2<=1/(b-2) give by Cauchy--Schwarz
\[
{\lambda\over P}\ge\sqrt{{b-2\over2P_0}}. \tag{13}
\]

The needed interquartile gap follows directly from a uniform local-limit
lemma.  An increment of `T=Gamma(b,1)+c Bin(b,p_c)` is iid with a uniform
Cramer gap for `c` in a compact neighbourhood of `mu`.  Fourier inversion
therefore gives, uniformly on every fixed standardized compact set,
\[
\left|\sqrt{bv(c)}f_{b,c}(bm(c)+\sqrt{bv(c)}x)-\varphi(x)\right|
\le C_Kb^{-1/2}. \tag{14}
\]
Near zero this is the cumulant expansion
`log psi(u/sqrt(b))=-u^2/2+O(|u|^3/sqrt(b))`; outside a fixed neighbourhood
the Cramer gap is exponentially small.  Chebyshev plus the positive upper
and lower central bounds in (14) imply
\[
\log f_{b,c}(T)+\tfrac12\log b=O_p(1),\quad
W-(T_1+T_2)+\log b=O_p(1),\quad
{W-Nm(c)+\log b\over\sqrt{Nv(c)}}\Rightarrow N(0,1). \tag{15}
\]
The last convergence is the triangular iid CLT from the same characteristic
function expansion.  It is uniform for `c_N->mu`, so the first and third
quartiles of `W` are separated by `c sqrt(N)`.

Since r054 gives `P_0=o(1)`, the KKT level `ell=-log lambda` lies below the
lower quartile.  On the upper quartile, the elementary lower profile bound
\(g(u)\ge\min\{(u-\log N)_+/N,1\}\) yields
\[
P_0=Eg(W-\ell)\ge cN^{-1/2}. \tag{16}
\]

Combining (12), (13), and (16), for all large `N`,
\[
 {\int q^-d\xi\over P}\ge c\sqrt{NP_0}\ge cN^{1/4},
\qquad
 {\int|q|d\xi\over P}=1+2{\int q^-d\xi\over P}\ge cN^{1/4}. \tag{17}
\]

Normalization and the mean-one Fenchel identity give integral q dxi=P, which justifies the total-variation equality above. This is the useful conclusion: the normalized signed density `q dξ/P` has
diverging total variation.  It cannot be treated as a probability law, and
the contrast moment in (6) cannot justify a positive-measure variance or a
local probability argument without a genuinely new sign/localization input.
It gives no new upper bound for a score and no general `C_N` statement.

## Computation record

`dilation55.wl` ran in Wolfram 15.0.0 for Mac OS X ARM.  It verified (1),
the parametric derivative identity, the 64-dimensional root, and the exact
moments/tail formula in (9).  The displayed consequences are analytic and
are conditional only on the cited r054 regular-KKT/profile estimates; the
uniform local-limit step used for (16) is stated and proved in (14)--(15).
