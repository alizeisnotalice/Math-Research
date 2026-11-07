# A self-contained local moderate-deviation lemma for Gamma--Bernoulli sums

This supplies the precise density/CDF input used in efficiency61. It does
not rely on the proof skeleton in r016. All constants below can be chosen
uniformly for `c` in any fixed compact subset of `(0,infinity)`.

Let `Z_c=E+cJ`, with `E~Exp(1)`, `J~Bernoulli(p_c)`,
`p_c=(1+exp(c))^(-1)`, independent. Write `K_c(theta)=log E exp(theta Z_c)`,
`m_c=K_c'(0)`, `sigma_c^2=K_c''(0)`, and let `F_{n,c}, f_{n,c}` be the
CDF and density of the sum of `n` independent copies. For fixed `0<a<b`,
uniformly for `a sqrt(log n)<=z<=b sqrt(log n)`, put
`x=n m_c-sigma_c sqrt(n) z`. Then
\[
 f_{n,c}(x)=\frac{e^{-z^2/2}}{\sigma_c\sqrt{2\pi n}}
 \{1+O((1+z^3)/\sqrt n)\},                                      \tag{1}
\]
\[
 F_{n,c}(x)=\frac{e^{-z^2/2}}{\sqrt{2\pi}z}
 \{1+O(z^{-2}+(1+z^3)/\sqrt n)\}.                               \tag{2}
\]
In particular, `f/F=(z/(sigma_c sqrt n))(1+o(1))` in this range.

## 1. A uniform local CLT under a small exponential tilt

The exact MGF is
\[
 e^{K_c(\theta)}=\frac{1-p_c+p_c e^{c\theta}}{1-\theta},\qquad\theta<1.
\]
Fix `theta0<1/2`. Under tilt `theta` with `|theta|<=theta0`, the law is
again a sum of an exponential of rate `1-theta` and an independent
Bernoulli multiple `c`, now with probability
`p_{c,theta}=(1+exp((1-theta)c))^(-1)`. Its variance is bounded above and
away from zero; its absolute moments through order three and a small
exponential moment are uniformly bounded on the stated compact set.

Let `chi_{c,theta}` be the centered characteristic function and
`sigma_{c,theta}^2=K_c''(theta)`. Its explicit exponential factor gives
\[
 |\chi_{c,\theta}(t)|\le
 \left(1+\frac{t^2}{(1-\theta)^2}\right)^{-1/2}.                 \tag{3}
\]
Also, uniformly for small real `t`,
`log chi(t)=-sigma_{c,theta}^2 t^2/2+O(|t|^3)` and
`|chi(t)|<=exp(-c0 t^2)`, with fixed positive `c0`.

These statements imply the following absolute local CLT for the tilted
sum, centered at its own mean:
\[
 \sup_y\left|\sigma_{c,\theta}\sqrt n\,
 f_{n,c,\theta}(nK_c'(\theta)+\sigma_{c,\theta}\sqrt n\,y)
 -\varphi(y)\right|\le C n^{-1/2}.                             \tag{4}
\]
Here are the integral estimates, so no local-limit theorem is being assumed.
Fourier inversion applies for `n>=2` by (3). For the normalized Fourier
variable `t`, split its integral at `n^(1/8)` and `epsilon sqrt n`.
On the first part, the logarithmic expansion and the small-argument
modulus bound imply
\[
 |\chi(t/(\sigma\sqrt n))^n-e^{-t^2/2}|
 \le C n^{-1/2}|t|^3 e^{-c_1t^2}.
\]
Its integral is `O(n^(-1/2))`. On the middle part both integrands have
integral `O(exp(-c2 n^(1/4)))`. On the remaining part, (3) bounds the
integral by
\[
 C\sqrt n\int_\epsilon^\infty(1+c_3u^2)^{-n/2}du
 \le C\sqrt n\,\rho^n,
 \qquad 0<\rho<1,
\]
for large `n`: retain half the exponent to get the factor `rho^n`, and
integrate the other half, bounded by an integrable fixed power. These bounds
are uniform in `c,theta` and prove (4), including a uniform constant.

## 2. Tilt at the target point

For `x` in the displayed moderate range, there is a unique small negative
`theta` solving `K_c'(theta)=x/n`, by strict convexity. Taylor expansion
with uniformly bounded third derivative gives
\[
 \theta=-\frac{z}{\sigma_c\sqrt n}+O(z^2/n),\quad
 \sigma_{c,\theta}=\sigma_c(1+O(z/\sqrt n)),
\]
\[
 I=x\theta/n-K_c(\theta)
 =\frac{z^2}{2n}+O(z^3/n^{3/2}).                                \tag{5}
\]
The exact change of measure at the mean of the tilted sum yields
\[
 f_{n,c}(x)=e^{-nI} f_{n,c,\theta}(x).
\]
Apply (4) at `y=0`; since `z^3/sqrt n=o(1)`, (5) proves (1).

For the CDF the same change of measure is exact:
\[
 F_{n,c}(x)=e^{-nI}\int_0^\infty e^{\theta v}
 f_{n,c,\theta}(x-v)dv.                                        \tag{6}
\]
Put `alpha=|theta|`. Substituting (4) in (6) has absolute error at most
`C/(n alpha)`. The Gaussian part equals
\[
 \frac{1}{\sigma_{c,\theta}\sqrt{2\pi n}\,\alpha}
 \left\{1+O\left((n\alpha^2)^{-1}\right)\right\},               \tag{7}
\]
because after `w=alpha v` its integral is
`int_0^infinity e^(-w) exp[-w^2/(2n sigma_{c,theta}^2 alpha^2)]dw`;
`0<=1-exp(-r)<=r` bounds its difference from one. The local-CLT error
relative to (7) is `O(n^(-1/2))`. Since
`alpha sigma_{c,theta} sqrt n=z(1+O(z/sqrt n))`, (5)--(7) prove (2).

## 3. The particular quantile and a fixed logarithmic neighborhood

Set `c=mu`, `m=N-2`, `A_N=1+sqrt(log N/N)`, and let
`F_{m,mu}(q_N)=1-A_N^(-1)`. Here `m_mu=sigma_mu^2=mu`.
At `z=(1/2)sqrt(log m)` and `z=2sqrt(log m)`, (2) brackets this tail
probability for all sufficiently large `N`. The actual quantile is
therefore inside that moderate range, without presupposing its location.
Taking logarithms in (2) gives
\[
 z(q_N)^2=\log N-2\log\log N+O(1),\qquad
 f_{m,\mu}(q_N)=\frac{\log N}{\sqrt\mu\,N}(1+o(1)).              \tag{8}
\]
For every fixed `C`, (1) on the fixed neighborhood
`|w|<=C log N` gives
\[
 \sup_{|w|\le C\log N}
 \left|\frac{f_{m,\mu}(q_N+w)}{f_{m,\mu}(q_N)}-1\right|=o(1).  \tag{9}
\]
Indeed its Gaussian log-density change is
`O((log N)^(3/2)/sqrt N)`, and the relative errors in (1) are uniformly
o(1). Equations (8)--(9) supply a completed, independent basis for the
noncircular CDF bracketing in efficiency61; they assert no endpoint theorem.
