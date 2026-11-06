# Uniform common-parameter curvature on the regular window

Scope: full positive quadrant, one receiver, two equal Gamma–Bernoulli blocks N=2b, cap two, target one. Fix K and assume |c−mu|<=K sqrt(logN/N). Derivatives are along (c1,c2)=(c,c). All differentiations/coarea interchanges satisfy the regular KKT hypotheses of r050/r054. This is not a theorem about the strong maximal constant C_N.

Let P be the minimized actual Phi_N budget, lambda its multiplier, H=h(lambda R), nu=R dxi, r=lambda/P, and L=logN. Use the complete nonnegative response form
\[
B(f,g)=\int_{\{\log(\lambda R)=0\}}fg\,d\nu/|\nabla\log(\lambda R)|
 +\int_{N<\lambda R<\Phi_N'(2)}{\lambda R\over\Phi_N''(H)}fg\,d\nu.
\]
In each block set s_i=partial_c logR_b(T_i), k_i=partial_t logR_b(T_i), s_+=s_1+s_2, k_+=k_1+k_2. Write B00=B(1,1), B0s=B(1,s_+), Bss=B(s_+,s_+), M=E_nu Hs_+ and T=E_nu H R_{++}/R.

We prove uniformly in this window
\[
\boxed{(\log P)''(c,c)=-Nm''(c)+O_K(N^{3/4}L^3)
=Np_\mu+O_K(N^{3/4}L^3).} \tag{1}
\]
No stationarity assumption is made. Thus on any connected C² regular window the log-budget is strictly convex for sufficiently large N. Existence of a stationary point is not inferred.

## 1. Improved sublevel mass and mean score

Set W=logR_b(T1)+logR_b(T2), ell=−loglambda and pi0=Pr_nu(W<ell). Since the KKT response is zero below the lower level and at least one above it, EH=1 gives
\[
\pi_0=E(H-1)_+=Eg(W-\ell),\quad E|H-1|=2\pi_0,\quad
\min\{(u-\log N)_+/N,1\}\le g(u)\le6(u-\log N)_+/N.
\tag{2}
\]
There are no level atoms for these block kernels.

Here are details of the concentration step. Let w0=Nm(c)−logb. On the moderate band |Ti−bm(c)|<=D sqrt(bL), compact exponential tilting and the Fourier density lower bound at the tilted mean give
\[
\log f_{b,c}(T_i)\ge-\tfrac12\log b-C_D L.
\]
The tilt parameter is O(sqrt(L/b)), the rate cost is O(L), and the tilted centered density is bounded below by a positive constant times b^(-1/2). These are exactly the compact-tilt Fourier estimates of r026/r053; the tilted exponential rate and Bernoulli parameter remain in fixed compact sets. Uniform Chernoff bounds make the two-block band complement at most N^(-A), for any fixed A, by increasing D.

Consequently W>=w0−C sqrt(NL) on that band. If ell were below w0−C sqrt(NL), then Pr(W<ell)<=N^(-A), hence by (2) Eg<=N^(-A). On a fixed central band of probability bounded below, the central density lower bound instead gives W>=w0−C' sqrtN−C'L. The lower estimate in (2) would make Eg at least a constant times sqrt(L/N), a contradiction for A>1. Thus
\[
\ell\ge w_0-C_K\sqrt{NL}.
\]
The global Gamma-mixture maximum f_b<=C/sqrtb gives W<=T1+T2−logb+C. The upper estimate in (2) and E|T1+T2−Nm|=O(sqrtN) now yield
\[
\pi_0=O_K(\sqrt{L/N}),\quad
\|H-1\|_2=O_K(N^{-1/4}L^{1/4}). \tag{3}
\]
On the same moderate band the r054 global/local score bound gives |s_+|<=C_K L. Outside it, |s_+|<=C_K N and the probability is arbitrarily polynomially small. Since Es_+=0, (2) proves
\[
|M|=|E(H-1)s_+|=O_K(N^{-1/2}L^{3/2}). \tag{4}
\]

## 2. Exact Ward normalization of B00

Location differentiation of the frozen source, including its lower jump, gives
\[
B(1,k_+)=E_\nu[\partial_1H+\partial_2H]
=E_\nu[H(2-k_+)]=2P/\lambda=2/r. \tag{5}
\]
The middle equality is full-quadrant integration by parts; the last is the two one-block Ward identities. This is the step unavailable from an upper bound on B00 alone.

Let l_i=k_i−1. The Gamma location-score conditional-expectation identity gives El_i=0 and El_i²<=1/(b−2). Push nu_i forward by Ui=logR_b(Ti). Its density q_b has norm <=C/sqrtb by r053. The lower response is integration on U1+U2=−loglambda. For any such level, disintegrating the square and using (l1+l2)²<=2l1²+2l2² gives
\[
B^{low}(k_+-2,k_+-2)\le4\|q_b\|_\infty E l_1^2=O(N^{-3/2}).
\]
This argument may use conditional second moments given Ui, so it does not assume an unproved inverse for the pushforward. For these kernels strict monotonicity also follows from k_i>=(b−1)/t>0. The soft factor lambdaR/Phi''<=6/N and independence give the additional O(N^(-2)) bound. Hence
\[
B(k_+-2,k_+-2)=O(N^{-3/2}). \tag{6}
\]
By (5), positive-form Cauchy, and (6),
\[
|B00-1/r|\le C\sqrt{B00}\,N^{-3/4}.
\]
The r053 multiplier bound r<=2N is uniform in this window. Put x=rB00. Then |x−1|<=C N^(-1/4)sqrtx. Solving this elementary quadratic in sqrtx proves
\[
\boxed{rB00=1+O(N^{-1/4}),\quad 1/B00-r=O(rN^{-1/4}).} \tag{7}
\]

## 3. The previously unresolved Schur/gradient cancellation

The exact r050 formula is
\[
(\log P)''=-r(T+Bss)+r{(M+B0s)^2\over B00}-r^2M^2.
\]
The stationary assumption is unnecessary once (7) is known. Expanding the last two terms gives
\[
rM^2(1/B00-r)+2rM B0s/B00+rB0s^2/B00. \tag{8}
\]
The r053 sum-score response estimate is Bss=O_K(L²/sqrtN). Positive-form Cauchy gives B0s²<=B00 Bss. Using (4), (7), and r<=2N, the three terms in (8) are respectively
\[
O_K(N^{3/4}L^3),\quad O_K(N^{3/4}L^{5/2}),\quad O_K(\sqrt N L^2).
\]
All are o(N). In particular, discarding the gradient square instead of comparing it with the Schur term would not establish this conclusion.

Finally r054's exact Ward residual decomposition gives
\[
T=Nm''(c)/r+E[HZ_+],\quad EZ_+=0,\quad\|Z_+\|_2=O_K(L^2).
\]
By (3), r|E[HZ_+]|<=C N^(3/4)L^(9/4). Together with rBss=O(sqrtN L²), (8) proves (1). Since m''(mu)=−p_mu and |c−mu|<=K sqrt(L/N), replacing −Nm''(c) by Np_mu costs only O_K(sqrt(NL)).

The result remains an auxiliary single-receiver theorem. It supplies no uniform multi-receiver Hessian, strict-target extension, general C_N upper bound, or new endpoint coefficient.
