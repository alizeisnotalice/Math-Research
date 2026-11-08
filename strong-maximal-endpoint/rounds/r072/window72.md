# Wide radial windows as a local-cost screen

## Verdict and scope

For the screening identity, fix a reference action `c` in the compact critical window and a fixed
cap `0<=F<=A`, the wide-window argument is valid. Define

\[
 j_W[\Phi_N(F)](u,X,\sigma)
 =\int_{W_N}e^{-s}[K_c^s\Phi_N(F_{u+s})](X,\sigma)\,ds. \tag{1}
\]

Its physical strong-`L^1` norm is exactly `|W_N|`. Thus a window of length
`O(\sqrt{N\log N})` yields screening exponent `1/2+o(1)`, below
`2\alpha=(9-\mu)/12\approx0.64346`. The r069 forbidden-band capture then
applies to a receiverwise windowed profile after a negligible mean repair.

The capture coefficient alpha below is used only for the optimized deterministic reference c_N=mu−sqrt(mu)/2 sqrt(log N/N); arbitrary critical c has its own scalar calibration. This is not a general upper bound: it leaves open the variable-action to
common-reference transfer and the common-anchor/source allocation.

## Choice of window

Let

\[
 p_{N,c}(s)={e^{-s}H_{N,c}(s)\over(1+e^{-c})^N}.
\]

Write `L=log N` and assume `|c−mu|<=D sqrt(L/N)` with D fixed. Each increment is `E+cB`, where E is Exp(1), B is Bernoulli with probability `p_c=1/(1+e^c)`, and they are independent. Its mean `m(c)=1+c p_c` satisfies `m(mu)=mu`, `m'(mu)=0`, and `|m''(c)|<=1` on `[0,2]`: use `|p'|<=1/4`, `|p''|<=1/4` in `m''=2p'+cp''`. Hence `|Nm(c)−Nmu|<=D²L/2`.

For `|t|<=1/2`, the centered exponential log moment is at most `t²`, by its power series. The centered Bernoulli log moment has second derivative at most `c²/4`, so it is at most `c²t²/8`. Thus for `c<=2`,

\[
 \log\mathbb E_c e^{t(S-Nm(c))}\le\tfrac32 Nt^2.
\]

For `0<d<=3N/2`, Chernoff with `|t|=d/(3N)` gives

\[
 \mathbb P_c(|S-Nm(c)|>d)\le2e^{-d^2/(6N)}.
\]

Take `W_N=[Nmu−7 sqrt(NL),Nmu+7 sqrt(NL)] intersect [0,infinity)`. For all sufficiently large N (in particular `N>=16L` and `D²L/2<=sqrt(NL)`), the distance from its endpoints to the mean is at least `6 sqrt(NL)`. Therefore

\[
 \boxed{\sup_c\mathbb P_c(S\notin W_N)\le2N^{-6}.}\tag{2}
\]

No Gaussian approximation or floating-tail certificate is used in this bound. More generally its constant can be enlarged for any prescribed polynomial tail power.

The constant 7 already ensures, for all large N at the calibrated references, that the forbidden band

\[
 I_N=[q_\circ-K\log N,\,q_\circ+K\log N]                 \tag{3}
\]

is contained in `W_N`. This is possible since
`N\mu-q_\circ=O(\sqrt{N\log N})` and `K\log N` is smaller.

## Exact physical `L^1` calculation

For every nonnegative `G`, Tonelli, physical radial translation, and joint
label Markovness give

\[
\begin{aligned}
 \|j_WG\|_{L^1(e^{-u}du\,dH\,d\lambda)}
 &=\int_{W_N}e^{-s}\int e^{-u}[K_c^sG_{u+s}](X,\sigma)
       \,du\,dH\,d\lambda\,ds\\
 &=\int_{W_N}\|G\|_{L^1(e^{-u}du\,dH\,d\lambda)}\,ds\\
 &=|W_N|\,\|G\|_1.                                      \tag{4}
\end{aligned}
\]

There is no label factor missing in (4): `K_c^s` is Markov on the full joint
label space. The identity fails for the full interval `(0,\infty)` only
because its length is infinite.

Set `G=\Phi_N(F)`, `\mathscr J=\|G\|_1`, and
`t=J_{*,2}N^\varepsilon`. Markov's inequality gives

\[
 \nu\{j_WG>t\}
 \le |W_N|{\mathscr J\over J_{*,2}N^\varepsilon}
 =Q_{\rm ref}N^{1/2-\varepsilon+o(1)}\mathscr J.          \tag{5}
\]

Choose

\[
 {1\over2}<\varepsilon<2\alpha,\qquad K>\alpha.          \tag{6}
\]

Then (5) is power-smaller than the reference scale and
`\varepsilon<2K` still holds as `K\downarrow\alpha`. Hence the capture
loss on the good set is `N^{K-\alpha+o(1)}` with arbitrarily small fixed
exponent.

## Receiverwise truncation and scalar capture

At a fixed receiver put

\[
 h(s)=[K_c^sF_{u+s}](X,\sigma),
 \qquad h_W(s)=h(s)\mathbf1_{W_N}(s).                      \tag{7}
\]

This does **not** alter the global source `F`: it is only a receiverwise
scalar profile, which the scalar capture theorem is free to use. If the
receiver is successful at the fixed reference action, then `\int p_{N,c}h
\ge1`. Thus (2) and `0\le h\le A` imply

\[
 \int p_{N,c}h_W\ge1-\eta_N,
 \qquad\eta_N=2A N^{-6}.                                   \tag{8}
\]

In particular `\eta_N=o(\log N/N)`. Rescale
`\widetilde h=h_W/(1-\eta_N)`. It has mean at least one, its cap is below a
slightly larger fixed cap, and

\[
 J(\widetilde h)\le(1-\eta_N)^{-N}J(h_W),
 \qquad\log(1-\eta_N)^{-N}=O(N^{-5})=o(1).               \tag{9}
\]

Conditional Jensen gives

\[
 J(h_W)=\int_{W_N}e^{-s}\Phi_N(h(s))ds
 \le j_W[\Phi_N(F)](u,X,\sigma).                          \tag{10}
\]

Consequently on a good receiver the r069 forbidden-band theorem applies to
`\widetilde h`. Since `I_N\subset W_N`, its captured cost is untruncated (now specialize to the optimized shifted reference). Explicitly, on I_N,

\[
 J_{I_N}(\widetilde h)\le(1-\eta_N)^{-N}j_{I_N}[\Phi_N(F)],
 \quad j_{I_N}[\Phi_N(F)]\ge(1-\eta_N)^N J_{I_N}(\widetilde h).
\]

The factor is 1−o(1) for the pure tail repair and N^{o(1)} for the reference-mean extension. Hence the scalar captured lower bound transfers to the original common-source band cost:

\[
 \int_{I_N}e^{-s}[K_c^s\Phi_N(F_{u+s})](X,\sigma)ds
 \ge J_{*,2}N^{\alpha-K+o(1)}.                             \tag{11}
\]

The repair in (9) is far below the permitted
`o(\log N/N)` mean error, and a fixed enlarged cap is already covered by
the r069 cap repair.

## Why full `j` is unnecessary here

The full cost `j=\int_0^\infty` has infinite physical `L^1` norm because
each radial slice costs one after the physical change of variables. Scalar
capture needs only a near-unit mean under `p_{N,c}` and an Orlicz cost bound
for the profile being captured. Equations (8)--(10) provide both for
`h_W`, while (3) ensures that capture on `I_N` is the actual physical band
cost. Thus full `j` is not needed in this fixed-reference, capped setting.

The same proof works when reference mean is at least 1−eta_ref,N with eta_ref,N=o(log N/N): use total repair eta_ref,N+2A N^−6, whose cost is N^o(1). The windowed cost is still paid by the original single G. The conclusion remains conditional on this fixed-reference mean comparison. Window
screening does not transfer variable-action successful receivers to that
reference, and does not solve the moving-anchor/common-source gap.

