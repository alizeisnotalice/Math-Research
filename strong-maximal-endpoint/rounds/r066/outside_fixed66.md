# Fixed far-action localization: repaired proof

## The scoped statement

Let \(r\) and \(A\) be fixed (enlarge the cap to \(\max\{A,2\}\) if needed), \(N=rb\), \(L=\log N\), and
\(\rho=\sqrt{L/N}\).  For arbitrary deterministic \(c_i\ge0\), let
\[
 T_i=\sum_{j=1}^b Z_{c_i,j},\qquad
 Z_c={\rm Exp}(1)+c\,{\rm Bernoulli}((1+e^c)^{-1}),
\]
with all summands independent.  If \(0\le H\le A\) and
\[
 J_+(H)\le e^{-s},\qquad s=N\mu-d_N+O(L),\qquad
 {\mathbb E}_{\mathbf c}H(T)\ge1,                         \tag{1}
\]
then
\[
 \max_{i\le r}|c_i-\mu|=O(\rho).                           \tag{2}
\]
The result is only for fixed equal blocks and a fixed cap.  It does not
make a claim about physical action-dependent angular kernels, arbitrary
rectangles, or the general strong maximal constant.

Here Φ_N(h)=h(1+log⁺h)^(N−1), J_+(H)=∫_{R_+^r}e^(−Σt_i)Φ_N(H(t))dt, and μ>1 is the unique root μ=1+e^(−μ). Define z_N>0 by barΦ(z_N)=√(μ/N)z_N and d_N=√(μN)z_N. Mills' inequalities give z_N~√(log N), hence d_N~√(μN log N). Applying (5) below at the common-μ threshold shifted by O(L) gives F_μ(Nμ−d_N+O(L))=d_N/N+O(L²/N); this also follows from r063/refined63.md.

Write
\[
 M=\sum_i b\,m(c_i),\quad V=\sum_i b\,\nu(c_i),\quad
 \Delta=N\mu-M,\quad
 a^2=\frac1r\sum_i(c_i-\mu)^2,                             \tag{3}
\]
where
\[
 m(c)=1+\frac c{1+e^c},\qquad
 \nu(c)=1+\frac{c^2e^c}{(1+e^c)^2}.
\]
The elementary identities \(m(\mu)=\nu(\mu)=\mu\), together with the
unique global maximum of \(m\), give
\[
 \mu-m(c)\ge\kappa\min\{(c-\mu)^2,1\}\quad(c\ge0).          \tag{4}
\]

## 1. Uniform Cramér estimate, including its Fourier proof

We use the following lemma in the range in which its relative error becomes
the required \(O(L^2/N)\) *absolute* CDF error.

**Lemma.**  Let \(S=\sum_{j=1}^N Z_{c_j}\) for any
\((c_j)\subset[0,\infty)\), with mean \(M\) and variance \(V\).  Uniformly
when
\[
 x=\frac{M-q}{\sqrt V}\in[1,C\sqrt L],
\]
\[
 {\mathbb P}\{S\le q\}
 =\overline\Phi(x)\left(1+O\!\left(\frac{(1+x)^3}{\sqrt N}\right)\right).
                                                                    \tag{5}
\]

**Proof.**  Let \(K(\theta)=\log{\mathbb E}e^{\theta S}\).  For a fixed
\(\theta_0<1/2\), under any \(\theta\in[-\theta_0,0]\), with
\(\gamma=1-\theta\), a tilted summand is exactly
\[
 {\rm Exp}(\gamma)+c\,{\rm Bernoulli}((1+e^{\gamma c})^{-1}). \tag{6}
\]
Thus its centered moments of every fixed order are uniformly bounded and
its variance is bounded above and below by positive constants.  Its
centered characteristic function factors into an exponential and a
Bernoulli factor, with
\[
 |\chi_{c,\theta}(t)|
 \le\frac{\gamma}{\sqrt{\gamma^2+t^2}}.                    \tag{7}
\]
The bound is uniform in \(c\): the Bernoulli factor has modulus at most
one.

Let \(V_\theta=K''(\theta)\asymp N\), \(\sigma_\theta=\sqrt{V_\theta}\),
and \(f_\theta\) be the tilted density. Choose the negative tilt by \(K'(\theta)=q\); its existence and \(\theta=O(\sqrt{L/N})\) follow from the uniform positive variance and the assumed range of \(x\). For
\(|u|\le\varepsilon\sqrt N\), Taylor expansion of the logarithm of the
product characteristic function gives
\[
 \log\prod_{j\le N}\chi_{c_j,\theta}(u/\sigma_\theta)
 =-\frac{u^2}{2}+O(|u|^3/\sqrt N),                         \tag{8}
\]
and, after decreasing \(\varepsilon\), its real part is at most \(-cu^2\).
For the complementary frequencies, (7), raised to the \(N\)-th power,
gives the integrable bound \(\int_{|u|>\varepsilon\sqrt N}(1+u^2/(C N))^{-N/2}du\le C\sqrt N e^{-cN}\).  Fourier inversion therefore gives
\[
 \sup_y\left|\sigma_\theta f_\theta(q+\sigma_\theta y)-\varphi(y)\right|
 \le C N^{-1/2}.                                           \tag{9}
\]
Indeed, on the central range the Fourier-integral difference is at most
\[
 {C\over\sqrt N}\int_{\mathbb R}|u|^3e^{-cu^2}\,du,
\]
and the two complementary integrals just described are exponentially
small.  This is a uniform density local CLT, not an action-compactness
argument.

Choose the negative tilt satisfying \(K'(\theta)=q\), set
\[
 I=\theta q-K(\theta),\qquad\lambda=-\theta\sigma_\theta.
\]
The exact tilted Laplace identity is
\[
 {\mathbb P}\{S\le q\}
 =e^{-I}\int_0^\infty e^{\theta v}f_\theta(q-v)\,dv.       \tag{10}
\]
Using (9), the error in the integral is at most
\(C/(\lambda\sqrt N)\).  The normal integral is retained exactly:
\[
 \int_0^\infty e^{-\lambda t}\varphi(t)\,dt
 =e^{\lambda^2/2}\overline\Phi(\lambda)
 \asymp (1+\lambda)^{-1}.                                  \tag{11}
\]
Hence the relative Fourier error is \(O(N^{-1/2})\).  Uniform cumulant
Taylor expansion following from (6) gives
\[
 \lambda=x+O(x^2/\sqrt N),\qquad
 I-\lambda^2/2=O(x^3/\sqrt N).                             \tag{12}
\]
Since \(|(\log\overline\Phi)'(x)|=O(1+x)\), (10)--(12) imply (5).
\(\square\)

In particular, when \(x\asymp\sqrt L\), (5) has relative
error \(O(L^{3/2}/\sqrt N)\); since
\(\overline\Phi(x)\asymp\rho\), its absolute error is
\[
 O(\rho L^{3/2}/\sqrt N)=O(L^2/N).                         \tag{13}
\]

## 2. A quantitative budget envelope

Fix a sufficiently large constant \(K\) and put
\[
 q=s-KL.
\]
Let \(P=(H-1)_+\), \(D=(1-H)_+\), and define the clipped scalar envelope
\[
 U(v)=\min\{A-1,\,[\Phi_N^{-1}(e^{v_+})-1]_+\}.             \tag{14}
\]
Here \(v_+=\max(v,0)\).  The clipping is essential: an unclipped inverse
envelope need not have a usable expectation in every fixed block dimension.

Every block density is a convex mixture of translates of a Gamma density
of supremum at most one, so the joint density \(f_{\mathbf c}\le1\) for
all actions.  First,
\[
 {\mathbb E}[H\,1_{\{S\le q\}}]
 \le\int_{\{S\le q\}}H(t)\,dt
 \le e^qJ_+(H)\le N^{-K}.                                  \tag{15}
\]
Thus
\[
 ED\ge F_{\mathbf c}(q)-N^{-K}.                            \tag{16}
\]
On the set \(\{P>U(S-q)\}\), the cap in (14) cannot be active and
\(\Phi_N(H)>e^{S-q}\).  Its Lebesgue measure, and hence its
\({\mathbb P}_{\mathbf c}\)-measure, is at most
\[
 e^qJ_+(H)\le N^{-K}.                                      \tag{17}
\]
Consequently \(EP\le EU+A N^{-K}\).  Exact success, \(ED\le EP\), yields
\[
 F_{\mathbf c}(q)\le EU+(A+1)N^{-K}.                       \tag{18}
\]
This is the required quantitative substitute for an \(o(\rho)\) envelope
error.

For \(u\ge0\), the exact modular formula is
\[
 \log\Phi_N(1+u)
 =\log(1+u)+(N-1)\log(1+\log(1+u)).
\]
Its first derivative at \(0\) is \(N\), and its second derivative has
absolute value at most \(CN\) on a fixed interval \(0\le u\le u_0\).
Choose \(\delta>0\) so that the inverse of \(v=\log\Phi_N(1+u)\) stays
in that interval for \(0\le v\le\delta N\).  Taylor inversion then gives
\[
 U(v)\le \frac vN+C\frac{v+v^2}{N^2}.                     \tag{19}
\]
For \(v>\delta N\), \(U\le A-1\), while the uniform positive exponential
moment of \(Z_c\) and \(M\le N\mu\) give
\[
 {\mathbb P}\{S-q>\delta N\}\le e^{-cN}.                   \tag{20}
\]
Initially, put \(h=(M-q)_+\le d_N+O(L)\).  The uniform positive
exponential moment gives the usual centered sub-gamma bound
\[
 {\mathbb P}\{S-M>t\}\le
 C\exp\{-c\min(t^2/N,t)\}.
\]
Integrating it, and using
\((S-q)_+\le h+(S-M)_+\), gives
\[
 E(S-q)_+=O(\sqrt{NL}),\qquad E(S-q)_+^2=O(NL).             \tag{21}
\]
Equations (19)--(21) imply \(EU=O(\rho)\), and (18) gives
\[
 F_{\mathbf c}(q)\le C\rho.                                \tag{22}
\]

## 3. Mean localization, then sharp localization

The uniform Berry--Esseen consequence of the proof of (5) and (22) imply
\[
 M-q\ge c\sqrt{NL}.
\]
Since \(M\le N\mu\) and \(q=N\mu-d_N+O(L)\),
\[
 \Delta=O(\sqrt{NL}),\qquad a=O((L/N)^{1/4})               \tag{23}
\]
by (4): first \(\sum_i\min\{(c_i-\mu)^2,1\}=o(1)\), so each action enters a fixed neighborhood of \(\mu\), where the minimum in (4) equals the square.

In this intermediate range, \(V=N\mu+O(Na)\).  Let \(F_\mu(q)\) be the CDF
with every action equal to \(\mu\).  Applying (5) to both arrays and
using the normal-tail logarithmic derivative gives
\[
 \log\frac{F_{\mathbf c}(q)}{F_\mu(q)}
 \ge c\sqrt{\frac LN}\Delta-CLa-C\frac{L^{3/2}}{\sqrt N}. \tag{24}
\]
The first term is the mean loss, the second is the possible variance
change, and the last is the uniform Cramér remainder.  Also
\[
 F_\mu(q)=\frac{d_N}{N}+O(L^2/N)\asymp\rho.                \tag{25}
\]
By (22), \(La=o(1)\) in (23), so (24) forces
\[
 \Delta=O(\sqrt{N/L}),\qquad a=O((NL)^{-1/4}).             \tag{26}
\]

The same normal comparison now has absolute precision.  The normal CDF
derivative in the mean-loss direction is at least \(cL/N\), while its
variance displacement is at most \(C\rho La\).  Combining these facts
with (13) and (25) gives
\[
 F_{\mathbf c}(q)\ge\frac{d_N}{N}
  +c\frac LN\Delta-C\rho La-C\frac{L^2}{N}.                \tag{27}
\]

It remains to sharpen \(EU\).  From (26), the threshold \(q\) is below the
mean by a fixed multiple of \(d_N\).  Integrating (5) at \(q-y\) gives
\[
 F_{\mathbf c}(q-y)\le C\rho e^{-c\rho y}
 \quad(0\le y\le C\sqrt{NL});
\]
the remainder of the integral is exponentially bounded by the uniform
Chernoff estimate.  Hence
\[
 E(q-S)_+=\int_0^\infty F_{\mathbf c}(q-y)\,dy=O(1).       \tag{28}
\]
Therefore
\[
 E(S-q)_+=M-q+E(q-S)_+=d_N-\Delta+O(L),
\]
and the second-moment bound in (21) still holds.  Using (19)--(20) gives
\[
 EU\le\frac{d_N-\Delta}{N}+O(L/N).                         \tag{29}
\]
Now (18) improves to
\[
 F_{\mathbf c}(q)\le\frac{d_N-\Delta}{N}+O(L^2/N).         \tag{30}
\]
Combining (27), (30), and \(\Delta\ge\kappa Na^2\) gives
\[
 a^2\le C\rho a+C\rho^2,
\]
which proves (2).

## 4. Correct global sequencing

The \(R=C\sqrt{NL},\,W=N\) grid of global_blocks65.md presupposes
near-\(\mu\) actions.  For arbitrary actions,
\[
 E(T_i-T_r)=b\{m(c_i)-m(c_r)\}=O(N),
\]
so the original cut is invalid before (2).

Before localization one must instead take \(R=CN\) and \(W=N^2\).  The
same uniform positive exponential moment used in (20) gives
\[
 \sup_{\mathbf c\ge0}{\mathbb P}\{\|x(T)\|_\infty>CN\}
 \le e^{-\gamma N};
\]
the random-shift overlap is \(1+O(N^{-1})\).  After cell normalization,
the r026 point-budget estimate gives the lower radial cutoff and the
budget-weighted positive-tail argument gives an upper \(O(L)\) radial
cutoff (with the larger polynomial \(W^{r-1}\) absorbed by its fixed
constant).  On that strip the translated terminal satisfies (1), so the
result above first localizes its successful arbitrary action vector to a
fixed near-\(\mu\) box.  Only then may the r065 universal-action repair
and fixed-action block transfer be applied.

The resulting global corollary is explicit. For every bounded H:R^r→[0,A] with finite J_global(H)=∫e^(−Σw_i)Φ_N(H(w))dw,

\[
 \nu\{u:\sup_{\mathbf c\in[0,\infty)^r}E_{\mathbf c}H(u+T)>1\}
 \le(1+o(1))Q_NJ_{\rm global}(H),\qquad d\nu=e^{-\sum u_i}du.
\]

The o(1) is uniform over these H for fixed r and A. The broad localization, cellwise normalization, discarded budget-weighted radial tail, local action confinement, and r065 common-action repair must be applied in precisely that order. No part of this corollary removes within-block physical angular dependence or the fixed-cap/fixed-block restrictions.
