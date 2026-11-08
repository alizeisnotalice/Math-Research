# Conditional first action response: a dimension-free theorem

Let
\[
q_c(z)={e^{-z}(1+{\bf1}_{\{z>c\}})\over1+e^{-c}}{\bf1}_{\{z\ge0\}},
\qquad
\Phi_{b,c}(s;\xi)=
 {\mathbb E}_c\!\left[e^{i\sum_j\xi_jZ_j}\mid\sum_jZ_j=s\right].
\]
Fix \(K\), let
\[
d_b=\sqrt{\mu b}\,z_b,\qquad
\overline\Phi(z_b)=\sqrt{\mu/b}\,z_b,\qquad
|s-(b\mu-d_b)|\le C_0\log b,
\]
where \(C_0\) is fixed, and assume
\[
|c-\mu|\le K\sqrt{\log b/b},\qquad \sum_j\xi_j=0.
\]
Then, for all sufficiently large \(b\),
\[
\boxed{\quad
\sup_{\sum_j\xi_j=0}
|\partial_c\Phi_{b,c}(s;\xi)|\le C_K.\quad} \tag{T}
\]
All operator derivatives below are strong derivatives on each fixed L² input, represented by differentiated Fourier multipliers; the derivative operators have the displayed uniform norms.

Consequently the conditional angular convolution multiplier has a
dimension-free first \(c\)-derivative on \(L^2\) of the sum-zero
hyperplane. This is a conditional result; it does not invoke the
unconditioned multiplier bound.

More generally, let \(I,J\Subset(0,\infty)\) be fixed compact intervals.
For all sufficiently large \(b\), the conclusions (T), (17), and (23)
hold uniformly for
\[
c\in I,\qquad s/b\in J,                                  \tag{CP}
\]
with a constant \(C_{I,J}\).  The critical window above is the
specialization needed later; no extension to \(c\downarrow0\) or
\(c\to\infty\) is claimed.  In the compact statement, every later
constant denoted \(C_K\) is read as \(C_{I,J}\).

## 1. Exact normalized flux recursion

If \(G_{b,c}(s;\xi)\) is the simplex Fourier integral before division by
\(f_{b,c}(s)\), then, with \(p_c=(1+e^c)^{-1}\),
\[
\partial_cq_c=p_cq_c-p_c\delta_c.
\]
Thus
\[
\partial_c\Phi_{b,c}(s;\xi)=A_{b,c}(s)
\left[\Phi_{b,c}(s;\xi)-{1\over b}\sum_i
e^{i\xi_ic}\Phi_{b-1,c}(s-c;\xi_{\widehat i})\right],
\quad
A_{b,c}=bp_c{f_{b-1,c}(s-c)\over f_{b,c}(s)}. \tag{1}
\]
When \(\sum\xi_i=0\), the exact face recentering is
\[
e^{i\xi_ic}\Phi_{b-1,c}(s-c;\xi_{\widehat i})
=e^{i\alpha\xi_i}\Phi_{b-1,c}(s-c;\eta^{(i)}),
\quad
\eta^{(i)}_j=\xi_j+{\xi_i\over b-1},\quad
\alpha={bc-s\over b-1}. \tag{2}
\]
This identity is useful bookkeeping, but the proof of (T) below does not
need to assert a separate face-to-bulk estimate.

For the radial derivative, the labelled one-factor distributional identity
is
\[
\partial_zq_{c,\varepsilon}
=-q_{c,\varepsilon}+a_0\delta_0+\varepsilon p_c\delta_c,
\qquad a_0=(1+e^{-c})^{-1}. \tag{RF1}
\]
Define
\[
\begin{aligned}
{\rm Face}_0&={1\over b}\sum_i
\Phi_{b-1,c}^{\varepsilon_{\widehat i}}(s;\xi_{\widehat i}),\\
{\rm Face}_c&={1\over b}\sum_i\varepsilon_i e^{i\xi_ic}
\Phi_{b-1,c}^{\varepsilon_{\widehat i}}(s-c;\xi_{\widehat i}),\\
\alpha_0&=a_0{f_{b-1,c}(s)\over f_{b,c}(s)},\qquad
\alpha_c=p_c{f_{b-1,c}(s-c)\over f_{b,c}(s)}={A_{b,c}(s)\over b}.
\end{aligned}
\]
Applying (3) to one factor and quotienting gives the exact normalized
radial flux formula
\[
\partial_s\Phi_{b,c}^{\varepsilon}
=\alpha_0({\rm Face}_0-\Phi_{b,c}^{\varepsilon})
+\alpha_c({\rm Face}_c-\Phi_{b,c}^{\varepsilon}). \tag{RF2}
\]
For comparison, the labelled \(c\)-flux is exactly
\[
\partial_c\Phi_{b,c}^{\varepsilon}
=A_{b,c}(s)({\Phi}_{b,c}^{\varepsilon}-{\rm Face}_c). \tag{RF3}
\]
If \(K\) is the posterior number of jumps under the all-plus conditional
mixture and \(G=s-cK\), differentiation of the Gamma density gives the
check
\[
\alpha_0=(b-1){\mathbb E}[G^{-1}\mid S=s]-\alpha_c. \tag{RF4}
\]

## 2. A uniform shifted-product integral

Put
\[
w(t)=(1+t^2)^{-1/2},\qquad h(t)=\min(t^2,|t|).
\]
For every \(b\ge3\) and every real shifts \(x_1,\ldots,x_b\),
\[
\int_{\mathbb R}\left(\sum_i h(\tau+x_i)\right)
\prod_jw(\tau+x_j)\,d\tau\le {C\over\sqrt b}. \tag{3}
\]

First,
\[
h(t)w(t)\le(1+\sqrt2)(1-w(t)). \tag{4}
\]
For \(|t|\le1\), divide the two sides by \(w(t)\) and use
\(t^2/(\sqrt{1+t^2}-1)=\sqrt{1+t^2}+1\); the case \(|t|\ge1\) is
the same elementary rationalization.

Let independent nonnegative radii \(R_i\) have survival function
\(\Pr(R_i\ge r)=w(r)\), and write
\(I_i=[-x_i-R_i,-x_i+R_i]\). Layer-cake and (4) give
\[
\begin{aligned}
&\int\sum_i(1-w(\tau+x_i))\prod_{j\ne i}w(\tau+x_j)\,d\tau\\
&\hspace{8mm}=
{\mathbb E}\sum_i\bigl(|\cap_{j\ne i}I_j|-|\cap_jI_j|\bigr). \tag{5}
\end{aligned}
\]
For deterministic intervals, deleting an interval can enlarge their
intersection only if that interval supplies the unique maximal left
endpoint or the unique minimal right endpoint. Hence at most two summands
in (5) are nonzero. Each leave-one-out intersection has length at most
twice the second-smallest radius \(R_{(2)}\). Therefore the right side of
(5) is at most \(4{\mathbb E}R_{(2)}\).

The latter expectation is uniformly \(O(b^{-1/2})\), since
\[
{\mathbb E}R_{(2)}
=\int_0^\infty\{w(r)^b+b(1-w(r))w(r)^{b-1}\}\,dr. \tag{6}
\]
For \(0\le r\le1\),
\[
w(r)\le e^{-r^2/4},\qquad1-w(r)\le r^2/2,
\]
so the integral in (6) over this range is \(O(b^{-1/2})\). For \(r\ge1\),
it is bounded by
\[
(b+1)2^{-(b-3)/2}\int_1^\infty w(r)^2\,dr=O(b^{-1/2}).
\]
Equations (4)--(6) prove (3), uniformly in all shifts.

## 3. Mean-matched fixed tilt and drift cancellation

For the proof at a fixed \((b,c,s)\), choose the unique exponential tilt
\(\theta<1\), or equivalently rate \(\gamma=1-\theta>0\), such that the
tilted mean is \(m_{c,\theta}=s/b\).  Indeed,
\[
{s\over b}={m(\gamma c)\over\gamma}
={1\over\gamma}+{c\over1+e^{\gamma c}},
\qquad m(d)=1+{d\over1+e^d}.                              \tag{CP1}
\]
The right side is strictly decreasing in \(\gamma\).  Since
\(1\le m(d)\le\mu\), the compact condition (CP) gives
\[
{1\over\max J}\le\gamma\le{\mu\over\min J}.               \tag{CP2}
\]
Thus every tilted parameter used below remains compact.  In the critical
window, \(\theta<0\) and
\(\theta=O_K(\sqrt{\log b/b})\).  Keep this \(\theta\) fixed while
differentiating in \(c\). This is legitimate because conditioning on the
sum is invariant under every fixed exponential tilt; no derivative of
\(\theta\) is taken.

The tilted CF is
\[
\psi_{c,\theta}(t)
={\gamma\over\gamma-it}\,
 {1+e^{-\gamma c+ict}\over1+e^{-\gamma c}},
\qquad \gamma=1-\theta.
\]
Let \(m=m_{c,\theta}\), \(v=v_{c,\theta}\), and set
\(\phi(t)=e^{-imt}\psi_{c,\theta}(t)\). Uniformly in the present
parameter window,
\[
\partial_c\log\psi_{c,\theta}(t)=i(\partial_cm)t+R_c(t),
\qquad |R_c(t)|\le C_Kh(t), \tag{7}
\]
and
\[
\partial_t\log\psi_{c,\theta}(t)-im=-vt+R_t(t),
\qquad |R_t(t)|\le C_Kh(t),\qquad v\asymp1. \tag{8}
\]
These follow by Taylor expansion at zero; for \(|t|>1\), the tilted
Bernoulli denominator is uniformly separated from zero and the direct
global bounds are \(O_K(|t|)\).

For \(\sum\xi_i=0\), write
\[
P_\xi(\tau)=\prod_i\phi(\xi_i+\tau),\qquad
I_\xi=\int_{\mathbb R}P_\xi(\tau)\,d\tau,\qquad I_0=\int\phi(\tau)^b\,d\tau.
\]
At the fixed mean-matched point, the conditional CF is \(I_\xi/I_0\).
For its \(c\)-derivative, differentiate instead the equivalent *raw*
fixed-\(\theta\) Fourier numerator and denominator, then evaluate at this
point.  Denote these derivatives by \(\dot I_\xi,\dot I_0\).  This
distinction is needed because the centering factor \(e^{-imt}\) is used
only to represent the integrands at the evaluation point.  Since \(b\ge3\),
(3) and the tilted exponential factor justify integration by parts. From (8),
\[
0=\int P_\xi'(\tau)\,d\tau
=\int P_\xi(\tau)\left[-vb\tau+\sum_iR_t(\xi_i+\tau)\right]d\tau,
\]
so
\[
b\int\tau P_\xi(\tau)\,d\tau
={1\over v}\int P_\xi(\tau)\sum_iR_t(\xi_i+\tau)\,d\tau. \tag{9}
\]
The raw fixed-tilt \(c\)-score in the numerator is
\[
\sum_i\partial_c\log\psi_{c,\theta}(\xi_i+\tau)
=ib(\partial_cm)\tau+\sum_iR_c(\xi_i+\tau).
\]
Thus
\[
\dot I_\xi=\int P_\xi(\tau)
\sum_i\partial_c\log\psi_{c,\theta}(\xi_i+\tau)\,d\tau.
\]
Using (9), (7), (8), and (3), this is bounded by
\[
|\dot I_\xi|\le {C_K\over\sqrt b}. \tag{10}
\]
The same bound holds for \(\dot I_0\).

Finally, the uniform tilted local density CLT at its mean gives
\[
I_0=2\pi f^{(\theta)}_{b,c}(bm)
\asymp_K b^{-1/2}. \tag{11}
\]
For completeness, put \(u=\sqrt b\,t\) in
\(\int\phi(t)^b\,dt\). The uniform third-moment Taylor expansion gives the
Gaussian integral on bounded \(u\), while
\(\lvert\phi(t)\rvert\le\gamma/\sqrt{\gamma^2+t^2}\) makes the complement
integrable and negligible uniformly in the present tilt window.
This is the same uniform tilted-density Fourier argument recorded in
r066, outside_fixed66.md, Section 1; its hypotheses hold here because
\(\theta=O_K(\sqrt{\log b/b})\) and \(c\) stays in a fixed compact
positive interval.
Also \(|I_\xi/I_0|\le1\), because it is a conditional characteristic
function. Differentiating the raw fixed-tilt ratio and applying
(10)--(11) proves (T).

## 4. Joint sign labels

For \(\varepsilon\in\{+1,-1\}\), put
\[
q_{c,\varepsilon}(z)=
{e^{-z}(1+\varepsilon{\bf1}_{\{z>c\}})\over1+e^{-c}}{\bf1}_{\{z\ge0\}}.
\]
Then
\[
\partial_cq_{c,\varepsilon}
=p_cq_{c,\varepsilon}-\varepsilon p_c\delta_c,\qquad
\psi_{\varepsilon,c,\theta}(t)
={\gamma\over\gamma-it}{1+\varepsilon r e^{ict}\over1+r},
\quad r=e^{-\gamma c}. \tag{12}
\]
For a label vector \(\varepsilon=(\varepsilon_1,\ldots,\varepsilon_b)\),
\(\Phi^\varepsilon_{b,c}(s;\xi)\) denotes the corresponding product
transform divided by the all-plus conditional normalizer. The denominator
therefore remains the all-plus normalizer. Center every labelled
factor by the all-plus mean \(m_+=s/b\).  The plus residuals satisfy
(7)--(8), whereas direct differentiation gives, for the minus residuals
after subtraction of that same common drift,
\[
|R_{c,-}(t)|+|R_{t,-}(t)|\le C_K(1+|t|). \tag{13}
\]

Choose \(\delta=\pi/(2c_{\max})\), with the action window contained in
\((0,c_{\max}]\).  If \(|t|\le\delta\), then
\[
{|1-r e^{ict}|\over1+r}\le{\sqrt{1+r^2}\over1+r}\le\kappa<1 \tag{14}
\]
uniformly in the action and tilt window.  Call a minus coordinate with
\(|t|\le\delta\) marked-small.  At each integration point let \(m\) be
their number. Since \(m\kappa^m\) is uniformly bounded, their score
contribution is bounded by
\[
C\int_{\mathbb R}m(\tau)\kappa^{m(\tau)}
\prod_jw_\gamma(\tau+x_j)\,d\tau
\le C_Kb^{-1/2},\qquad
w_\gamma(t)={\gamma\over\sqrt{\gamma^2+t^2}}, \tag{15}
\]
because one first takes
\(\sup_{m\ge0}m\kappa^m<\infty\), and then Hölder bounds the remaining
shifted product integral by \(\|w_\gamma\|_b^b=O_K(b^{-1/2})\).

For a minus coordinate outside that small interval,
\[
(1+|t|)w_\gamma(t)\le C_{K,\delta}(1-w_\gamma(t)). \tag{16}
\]
The rescaled random-interval argument proving (3) bounds the sum of
marked-large terms. It also bounds all unmarked plus terms through their
\(h\)-residuals; dropping the remaining factor \(\kappa^m\le1\) only
enlarges the estimate. Therefore the labelled product obeys the analogue
of (9)--(10), with the common drift \(-v_+b\tau\).

For every \(\varepsilon\in\{\pm1\}^b\),
\[
\boxed{\quad
\sup_{\sum_j\xi_j=0}
|\partial_c\Phi^\varepsilon_{b,c}(s;\xi)|\le C_K.\quad} \tag{17}
\]
The quotient has modulus at most one: under the all-plus conditional law,
its numerator is the expectation of the Fourier character times the
product of \(\pm1\) jump labels.  Plancherel therefore gives the same
dimension-free first-response \(L^2\) multiplier bound uniformly over
joint sign characters.

## 5. Labelled operator and finite-radial comparison

Let \(H_b=\{x\in\mathbb R^b:\sum_i x_i=0\}\), and let
\(\mathcal K_c^s\) be the conditional angular convolution operator on
\(L^2(H_b\times\{\pm1\}^b)\).  Fourier transform on \(H_b\) and the
finite Fourier transform on the label group diagonalize it into the
multipliers \(\Phi^\varepsilon_{b,c}(s;\xi)\).  Hence (17) gives
\[
\|\partial_c\mathcal K_c^s\|_{2\to2}\le C_K. \tag{18}
\]
For every compact action interval \(J\) in the stated window, strong
differentiation, the fundamental theorem of calculus, and Minkowski give
\[
\left\|\sup_{c\in J}
\left|\mathcal K_c^sF-\mathcal K_{c_0}^sF\right|\right\|_2
\le C_K|J|\,\|F\|_2,\qquad c_0\in J. \tag{19}
\]
For fixed \(b,s\), the conditional joint kernel is bounded on its compact
simplex and is \(L^2\)-continuous in \(c\): its normalizer is positive,
and the moving jump hyperplanes have simplex measure zero, so dominated
convergence applies.  Convolution therefore has a natural jointly
continuous \((X,c)\) representative for every \(F\in L^2\).  One may first
take the supremum in (19) over a countable dense subset of \(J\), then use
this representative to identify it with the full supremum.  In
particular, the bound applies without a null-version ambiguity to an
arbitrary measurable policy \(c=c(X,\sigma)\in J\).

There is an exact fixed-radial physical-coordinate comparison.  Put
\(d=b-1\), take the barycentric receiver angle
\(\xi_0=b^{-1}{\bf1}\), and use
\[
X=V(\xi-\xi_0),\qquad Y=(V+s)(\eta-\xi_0).
\]
The affine conditional transport becomes
\[
Y=X+\Pi z,\qquad \Pi=I-b^{-1}{\bf1}{\bf1}^{\mathsf T}. \tag{20}
\]
Zero-extend the source terminal from its \(Y\)-simplex to \(H_b\).
On the receiver simplex \(S_V\), the conditional convolution of the
zero extension is exactly the physical angular average.  If \(B\) is
common to all source angles and labels, then it is fixed by every
conditional Markov kernel on \(S_V\); thus (18) applies to the
zero-extended \(F-B\).  The change of variables and (18) give the
finite-depth response estimate
\[
\|\partial_c\mathcal K_c^{V,s}F\|_{L^2(S_V,d\xi;\ell^2_\sigma)}
\le C_K\left({V+s\over V}\right)^{d/2}
\|F-B\|_{L^2(S_{V+s},d\eta;\ell^2_\sigma)}. \tag{21}
\]
After squaring, the finite-radial source/receiver Jacobians cancel exactly:
\[
V^d\|\partial_c\mathcal K_c^{V,s}F\|_{L^2(d\xi;\ell^2_\sigma)}^2
\le C_K^2(V+s)^d\|F-B\|_{L^2(d\eta;\ell^2_\sigma)}^2. \tag{22}
\]
For an action interval \(I\) in the same compact parameter window,
\[
\left\|\sup_{c\in I}
\left|\mathcal K_c^{V,s}F-\mathcal K_{c_0}^{V,s}F\right|\right\|_
{L^2(S_V,d\xi;\ell^2_\sigma)}
\le C_K|I|\left({V+s\over V}\right)^{d/2}
\|F-B\|_{L^2(S_{V+s},d\eta;\ell^2_\sigma)}. \tag{22a}
\]
This is a fixed-\((V,s)\) \(L^2\) interface, not a global source-profile
or receiver-cone allocation theorem.

## 6. Conditional radial sensitivity

At fixed \(c\), the same argument gives a stronger derivative in the
conditioned radial sum:
\[
\boxed{\quad
\|\partial_s\mathcal K_c^s\|_{L^2(H_b\times\{\pm1\}^b)\to
L^2(H_b\times\{\pm1\}^b)}
\le {C_K\over b}.\quad} \tag{23}
\]
This holds while \(s\) remains in the preceding critical window, and
includes all joint sign characters.

Indeed, use the same fixed mean-matched tilt only as a representation at
the evaluation point.  The derivative of the raw Fourier numerator with
respect to its conditioning value is
\[
\partial_sN_\xi(s)=-i\int_{\mathbb R}\tau P_\xi(\tau)\,d\tau.
\]
Equation (9), the residual bound, and the shifted-product lemma give
\[
\left|\int\tau P_\xi(\tau)\,d\tau\right|
\le {C_K\over b^{3/2}}. \tag{24}
\]
The same estimate holds for the all-plus denominator and for every
labelled numerator, using the marked-small/marked-large bounds in Section
4.  Dividing by \(I_0\asymp b^{-1/2}\) and differentiating the quotient
proves (23).

Thus, for any interval of critical values joining \(s\) and \(s+\delta\),
\[
\|\mathcal K_c^{s+\delta}-\mathcal K_c^s\|_{2\to2}
\le {C_K|\delta|\over b}. \tag{25}
\]
For a fixed angular source family \(\mathcal F(w,\cdot)\), this gives the
fibrewise radial-coherence comparison
\[
\begin{aligned}
&\|\mathcal K_c^{s+\delta}\mathcal F(u+s+\delta,\cdot)
-\mathcal K_c^s\mathcal F(u+s+\delta,\cdot)\|_2\\
&\hspace{45mm}\le {C_K|\delta|\over b}
\|\mathcal F(u+s+\delta,\cdot)\|_2. \tag{26}
\end{aligned}
\]
It compares only the two conditional kernels applied to the same
translated source fibre.  Controlling the separate source translation
\(\mathcal F(u+s+\delta)-\mathcal F(u+s)\), integrating (26), and charging
it to an actual source budget remain additional tasks.  Subject to such
control, a radial repair of size \(\delta\log N/N\) dominates this
kernel-only error \(\delta/N\) by a logarithmic factor for fixed blocks.

## 7. Singleton obstruction and remaining scope

The joint-label theorem deliberately requires \(b\ge3\). At fixed \(s>0\)
and frequency zero, a singleton block has
\[
\Phi^-_{1,c}(s;0)=
\begin{cases}
0,&c<s,\\
1,&c>s,
\end{cases} \tag{27}
\]
because below the observed value the two latent jump components have equal
weight, whereas above it only the no-jump component is possible. Thus no
uniform labelled \(C^1\) statement can include singleton blocks.

The theorem controls conditional Fourier response, including sign labels.
It does not by itself construct a compatible global source profile, carry
the fixed-\((V,s)\) comparison above into a global radial argument, or
allocate a single source deficit across many receiver cones. Those are
separate steps for a general strong-maximal upper bound.
