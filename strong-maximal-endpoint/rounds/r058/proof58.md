# r058: radial mixed coefficient

## Exact response representation

Scope is the full-positive-quadrant, frozen symmetric common-stationary target-one cap-two KKT source. Let N=2b, L=logN, r=lambda/P, and |c-mu|<=C_0 sqrt(L/N), with mu=1+exp(-mu) and p_mu=1/(1+exp(mu)). Constants depend only on the fixed window radius C_0. All r050/r054 regular KKT and integration assumptions are retained; statements are for sufficiently large N. Assume b>=6 for the local C4/implicit-function coefficient formula. Put \(V=\partial_{t_1}-\partial_{t_2}\), \(k_-=k_1-k_2\), \(s_+=s_1+s_2\), and \(s_{t,-}=\partial_{t_1}s_1-\partial_{t_2}s_2\). At the r057 base point,
\[
G_d=B(1,k_+)=2/r,\qquad G_{vv}=-B(k_-,k_-),\qquad G_{ad}=B(k_+,s_+).
\]
Here \(p_c(t_1,t_2)=f_{b,c}(t_1)f_{b,c}(t_2)\) is the probability density. Because \(V\log p_c=k_-\), one integration by parts gives
\[
G_{avv}=-B\bigl(k_-,s_{t,-}+s_+k_-\bigr).
\]
The lower-level term is understood distributionally: on smooth regions
\(VH_0=(\lambda R)h'(\lambda R)k_-\), while the jump is precisely the
lower-coarea part of the full response form \(B\).

With
\[
\beta={G_{ad}\over G_d}={B(k_+,s_+)\over B(1,k_+)},
\]
the radial correction has the exact centered form
\[
\boxed{K=-B\bigl(k_-,s_{t,-}+(s_+-\beta)k_-\bigr).} \tag{1}
\]
Thus \(K\) is not controlled by a Fisher bound on \(s\) alone: it also
requires the mixed location derivative \(s_{t,-}\).

## Curvature ratio

If \(G(\rho(a),0,a)=1\) and
\(C(a)=-G_{vv}(\rho(a),0,a)/G_d(\rho(a),0,a)\), stationarity gives
\[
C'(c)=-K/G_d,\qquad
\left[-{G_{vv}\over2G_d}\right]'_{a=c}=-K/(2G_d).
\]
So \(K\) is the common-action derivative of radially normalized angular
curvature.

## Established consequence from r057

The lower-coarea argument for \(l_+=k_+-2\) also applies to
\(k_-=l_1-l_2\). Using \(\|q_b\|_\infty=O_{C_0}(N^{-1/2})\),
\(El_i^2=O(N^{-1})\), and the soft factor \(O(1/N)\),
\[
B(k_-,k_-)=O_{C_0}(N^{-3/2}),\qquad
d_2=-G_{vv}/(2G_d)=O_{C_0}(N^{-1/2}).
\]

## A uniform \(K\) bound

On the moderate band \(\max_i|T_i-bm(c)|\le D\sqrt{bL}\), the r054 tilt
bridge gives
\[
|k_i-1|+|s_{i,t}|\le C_{C_0}\sqrt{L/N},\qquad |s_i|\le C_{C_0}L. \tag{7}
\]
The required transfer of this estimate to the response measure is available
for exactly the weights in (1). Write \(X=T-cJ\), \(m=b-1\),
\(\pi=(1+e^c)^{-1}\), and condition on \(T\). Then
\[
k=mE(X^{-1}\mid T),\qquad
{R_{ct}\over R}=b\pi mE(X^{-1}\mid T)-m(m-1)E(JX^{-2}\mid T). \tag{8}
\]
Together with the global \(|s|\le Cb\), this gives
\[
f_{b,c}|s_t|\le Cb\,(\hbox{mixtures of }\gamma_{b-1},\gamma_{b-2}). \tag{9}
\]
Let U=log R_b(T). Equation (8) gives k>0, so this pushforward has density q=f_{b,c}/k; functions in the following display are composed with its inverse. Since \(q=f_{b,c}/k\le tf_{b,c}/m\), \(t=X+cJ\), and
\(X\gamma_j=j\gamma_{j+1}\), for fixed \(j=0,1,2\)
\[
\sup_u q(u)k(u)^j=O_{C_0}(N^{-1/2}),\quad
\sup_u q(u)|s_t(u)|=O_{C_0}(\sqrt N),\quad
\sup_u q(u)k(u)|s_t(u)|=O_{C_0}(\sqrt N). \tag{10}
\]
The corresponding fixed moments of \(k,s_t\) are polynomial in \(N\), by
the same conditional inverse-Gamma formulas. Consequently, on a lower
coarea level, split according to whether either coordinate leaves the
moderate band. Use (10) on the other coordinate and Hölder--Chernoff on the
bad coordinate. Increasing \(D\) makes each of the finitely many weighted
tails appearing in (1) \(O_{C_0}(N^{-A})\), for arbitrary prescribed \(A\),
uniformly in the level. The soft part obeys the same conclusion using its
\(6/N\) factor and product moments. This is a finite-weight response-tail
lemma, not an assertion about all polynomial weights.

The Ward identity \(1/r=E[(H_0-1)(1-k_1)]\), (7), and
\(E|H_0-1|=O_{C_0}(\sqrt{L/N})\) now give
\[
r\gtrsim_{C_0} N/L,\qquad B_{00}\asymp r^{-1}=O_{C_0}(L/N). \tag{11}
\]
Initially \(|\beta|\le Cb\), because \(k_+>0\) makes it a response-weighted
average of \(s_+\). The tail lemma and (7) improve this to
\(\beta=O_{C_0}(L)\). Thus the integrand in (1) on the moderate band is
\[
\left|k_-\{s_{t,-}+(s_+-\beta)k_-\}\right|=O_{C_0}(L^2/N).
\]
Its response mass is \(B_{00}=O_{C_0}(L/N)\), while the weighted tails are
negligible. Hence
\[
\boxed{|K|=O_{C_0}(L^3/N^2).} \tag{12}
\]

Finally \(G_{aa}=T=Nm''(c)/r+E[H_0Z_+]\). The unscaled residual obeys \(E[H_0Z_+]=O_{C_0}(N^{-1/4}L^{9/4})=o(1)\) by r057 profile stability and its Z-plus bound. Thus \(r\le2N\) and \(m''(c)=-p_\mu+o(1)\) give
\(-G_{aa}\ge p_\mu/4\) for sufficiently large \(N\). Since \(G_d\ge1/N\),
the frozen quartic coefficient obeys
\[
\left|{K^2\over8G_{aa}G_d}\right|=O_{C_0}(L^6/N^3). \tag{13}
\]
This is only a coefficient bound. Uniform Taylor/implicit-function
remainder control is still required before using it on an angle interval
growing with \(N\).

## Symbolic audit

`analytic58.wl` was run with Wolfram Language 15.0.0 on 2026-11-07. For
generic positive \(p=e^{-t_1-t_2+r(t_1,a)+r(t_2,a)}\), all six algebraic
checks passed: the density jets, smooth \(VH\) response factor, one-IBP
product rule, centered response rearrangement, and curvature-ratio
derivative. This validates algebra only, not the lower-coarea tail lemma.


## Details behind the finite-weight tail lemma

The needed bounds use only the weights k^j (j<=2), |s_t| and k|s_t|, together with the global |s|<=Cb; no unrestricted polynomial-weight assertion is needed. If gamma_j is the unit-rate Gamma density, then gamma_b(x)/x=gamma_{b-1}(x)/(b-1) and gamma_b(x)/x²=gamma_{b-2}(x)/[(b-1)(b-2)]. Equation (8) and s_t=R_ct/R-sk therefore bound f_b|s_t| by Cb times normalized mixtures of shifted gamma_{b-1}, gamma_{b-2}. Their maxima are O(b^(-1/2)). Also t gamma_j(t-cJ)=[(t-cJ)+cJ]gamma_j, with (t-cJ)gamma_j=j gamma_{j+1} and J<=b. This proves both weighted |s_t| maxima in (10) after multiplication by t/(b-1) when necessary. For j=2, q k²=f_b k is exactly a mixture of gamma_{b-1}; j=1 is f_b and j=0 uses q<=t f_b/(b-1).

For every fixed p and sufficiently large b, conditional Jensen and E X^(-p)=Gamma(b-p)/Gamma(b) give ||k||_p=O_p(1), ||s_t||_p=O_p(b). Expand k_- s_{t,-} into four separable one-/two-block terms. The remaining terms in (1) are controlled by Cb times k_1²+k_2²+2k_1k_2, initially using |beta|<=Cb. On a lower coarea level, if block one is outside the moderate band, bound the other block's weighted pushforward density by (10); apply Holder to the first block's moment on the bad event. Interchange the roles for block two. The result is polynomial(b) times a positive power of the Chernoff tail, uniformly in the coarea level. Increasing the fixed band constant D makes this O(N^(-A)) for any prescribed A. The soft part follows from the same moments and its 6/N factor. These arguments also cover the weights k_+ s_+ used to bootstrap beta.

For clarity, the central s_t bound is a derivative of the r054 compact tilt bridge: s=bF(c,d)+O_{C^1}(1), F_d=−u'(d)/c with u'(mu)=0, and d_t=−d²/[bc v(d)]. Hence s_t=O(|d−mu|+1/b)=O(sqrt(L/N)). This is why a density-jet or first-mixed-score estimate, not a Fisher bound alone, is required.

The same concentration gives d_2=O_{C_0}(L/N): on the moderate band k_-²<=CL/N and B(k_-,k_-)<=CL B00/N plus a negligible tail, while rB00 stays bounded and d_2=r B(k_-,k_-)/4. Neither this nor the quartic coefficient bound controls the implicit-function remainder on an N-dependent angle interval.
