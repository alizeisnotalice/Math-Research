# r056: common-parameter curvature at a regular stationary point

Date: 2026-11-07. Scope: the r054 full-positive-quadrant, cap-two,
target-one, two-equal-block budget. This proves a uniform statement at
common-parameter stationary points, not throughout the parameter window.

## Exact common Hessian and its off-stationary obstruction

Let \(N=2b\), \(P(c)=P(c,c)\), and at a regular KKT point put
\[
 M_+=E_\nu[H(s_1+s_2)],\qquad
 T_{++}=E_\nu\left[H\left\{{R_{++}\over R}\right\}\right],
\]
where \(s_i=\partial_{c_i}\log R\).  Write
\(B_{0+}=B(1,s_1+s_2)\), \(B_{++}=B(s_1+s_2,s_1+s_2)\) for the complete
r050 coarea-plus-soft response form.  Its exact formula is
\[
 (\log P)''_+=-{\lambda\over P}(T_{++}+B_{++})
 \quad+{\lambda\over P}{(M_++B_{0+})^2\over B_{00}}
 -\left({\lambda M_+\over P}\right)^2. \tag{1}
\]
Also \(P'_+=-\lambda M_+\).  Away from stationarity, the last two terms
are a genuine Schur/gradient combination.  Existing r054 \(L^2\) estimates
do not bound it by \(o(N)\), so they do not justify a whole-window common
curvature theorem.

At a common stationary point \(M_+=0\).  Cauchy for the positive form \(B\)
gives \(B_{0+}^2\le B_{00}B_{++}\), and hence
\[
 (\log P)''_+=-{\lambda\over P}T_{++}
 -{\lambda\over P}\left(B_{++}-{B_{0+}^2\over B_{00}}\right). \tag{2}
\]

## Ward decomposition

Retain r054's one-block residual
\[
 e_i={R_{c_ic_i}\over R}-b m''(c)(1-k_i),\qquad
 Z_+=e_1+e_2+2s_1s_2. \tag{3}
\]
Then the common normalized second kernel derivative is
\[
 {R_{++}\over R}=b m''(c)(2-k_1-k_2)+Z_+. \tag{4}
\]
The r054 one-block Ward identity is
\(\lambda E_\nu H(1-k_i)=P\), therefore
\[
 T_{++}=Nm''(c){P\over\lambda}+E_\nu[H Z_+]. \tag{5}
\]
Equations (2) and (5) yield the exact stationary identity
\[
 \boxed{(\log P)''_+=-Nm''(c)
 -{\lambda\over P}E_\nu[H Z_+]
 -{\lambda\over P}\left(B_{++}-{B_{0+}^2\over B_{00}}\right).} \tag{6}
\]

Local Wolfram 15.0 checks archived in `code/common_checks56.wl` verified the plus/minus
decompositions in (4), together with \(m'(\mu)=0\) and
\(m''(\mu)=-p_\mu\).

## Uniform stationary curvature

Suppose \(|c-\mu|\le K\sqrt{\log N/N}\) and r054 regularity holds.  The
r054 estimates carry over from \(Z_-\) to \(Z_+\): block independence and
\(Es_i=0\) give \(EZ_+=0\), while changing the sign of \(2s_1s_2\) preserves
the relevant norm bounds. Thus
\[
 \|Z_+\|_2=O_K((\log N)^2),\qquad \|H-1\|_2=O(N^{-1/8}),
\]
and
\[
 |E_\nu[HZ_+]|=|E_\nu[(H-1)Z_+]|=
 O_K(N^{-1/8}(\log N)^2). \tag{7}
\]
The lower-coarea convolution and soft-response estimates used for
\(B_{--}\) apply with \(s_1+s_2\), whose second moment has the same order:
\[
 0\le B_{++}=O_K((\log N)^2/\sqrt N). \tag{8}
\]
As \(\lambda/P\le2N\), both error terms in (6) are \(o(N)\); and
\(-m''(c)=p_\mu+o(1)\) uniformly. Hence
\[
 \boxed{P'(c)=0\quad\Longrightarrow\quad
 (\log P)''(c)=Np_\mu+o_K(N).} \tag{9}
\]
All regular stationary points in the window are strict local minima of
\(\log P\). If \(P\) is \(C^2\) on a connected regular window, there is at
most one: every critical point crosses from negative to positive, and two
such crossings would force an intervening non-positive-curvature critical
point.

## Fixed-input implication for common receiver actions

Freeze the stationary KKT input \(H_0\), and define
\(A_{H_0}(c)=E_{(c,c)}H_0\).  Then
\[
 A_{H_0}'(c)=M_+=0,\qquad A_{H_0}''(c)=T_{++}.
\]
By (5), (7), and \(\lambda/P\le2N\),
\[
 {T_{++}\over NP/\lambda}=m''(c)+o_K(1)=-p_\mu+o_K(1). \tag{10}
\]
Since \(NP/\lambda\ge1/2\),
\[
 A_{H_0}''(c)\le-{p_\mu\over2}+o_K(1)<0. \tag{11}
\]
Thus a nearby receiver with common action \(c+\delta\), against the same
stationary optimal source, loses average at order \(\delta^2\).  At the
base receiver, varying its common action alone therefore has no
first-order gain and has a quadratic loss.  This does not rule out
angle-dependent actions at other receivers: mixed location/action terms
remain to be estimated.  This is not a claim for arbitrary fixed inputs.

## Local angular IFT corollary: the first possible adaptation gain is quartic

The preceding limitation has a precise local formulation.  Freeze the same
stationary symmetric input \(H_0\), and define
\[
 F(v,a)=E_{a,a}\,H_0(T_1+v,T_2-v). \tag{12}
\]
Here the two block increments have common parameter \(a\).  For \(b\ge6\),
the zero-extended Gamma--binomial density and bounded \(H_0\) give the
local \(C^4\) regularity needed below.  Exchange symmetry gives
\(F(v,a)=F(-v,a)\), while (11) gives
\[
 F_a(0,c)=0,\qquad F_{aa}(0,c)=T_{++}<0. \tag{13}
\]
The implicit-function theorem therefore supplies a unique local critical
action \(a(v)\), which is a local maximum and is even in \(v\).  In
particular
\[
 a'(0)=0,\qquad a''(0)=-{F_{avv}(0,c)\over F_{aa}(0,c)},
\]
and hence
\[
 a(v)=c-{F_{avv}(0,c)\over2F_{aa}(0,c)}v^2+o(v^2). \tag{14}
\]
Taylor expansion at fixed \(v\) then gives the exact leading adaptation
gain
\[
 F(v,a(v))-F(v,c)
 =-{F_{avv}(0,c)^2\over8F_{aa}(0,c)}v^4+o(v^4)\ge0. \tag{15}
\]
The displayed quartic coefficient is strictly positive precisely when the mixed derivative
\(F_{avv}(0,c)\) is nonzero.  This corollary neither asserts that it is
nonzero nor supplies a uniform \(v\)-window or an \(N\)-rate; computing
that mixed coefficient is the remaining common-action question.

## Limits

The target is exactly one. No strict-target extension, nonstationary common
curvature bound, or general strong-maximal conclusion is claimed.

The quartic comparison is between two averages at the same shifted receiver. That receiver may still have average below one. Strict acceptance, radial compensation, actual witness volume and the Orlicz denominator must be checked separately before any new endpoint lower bound is claimed.
