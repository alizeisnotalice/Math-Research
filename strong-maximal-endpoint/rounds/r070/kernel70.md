# The local-cost kernel: exact facts and the unresolved weak endpoint

Fix `N>=2`, an action `c>0`, and retain every Bernoulli label.  This note
uses Hausdorff measure on `H_N={x:\sum_i x_i=0}` and normalized Haar
probability on the label group.  It records exact operator facts only; it
does not assert a weak `(1,1)` theorem.

## Finite-action conditional kernel and physical local cost

Let

\[
 Z_i=E_i+cJ_i,\quad E_i\sim {\rm Exp}(1),\quad
 J_i\in\{0,1\},\quad P=I-N^{-1}{\bf1}{\bf1}^{\mathsf T},
\]

with Bernoulli probability P(J_i=1)=(1+exp(c))^(-1), and put `S=\sum_iZ_i`.
For an increment-label vector `J`, set

\[
 D_{s,J}=\{y\in H_N:y_i\ge cJ_i-s/N\},\qquad
 H_{N,c}(s)={1\over(N-1)!}\sum_{k=0}^N{N\choose k}(s-ck)_+^{N-1}.
\]

The joint conditional kernel of `(PZ,J)` given `S=s` has density, relative
to `dH` and counting measure in the increment label,

\[
 k_{c,s}(y,J)={\mathbf1_{D_{s,J}}(y)\over\sqrt N\,H_{N,c}(s)}. \tag{1}
\]

Thus the joint-label convolution `K_c^s` is Markov.  For example, with
`\chi_J` denoting the corresponding label-group increment, one can write

\[
 K_c^sG(X,\sigma)=\sum_J\int_{H_N}k_{c,s}(y,J)
 G(X+y,\sigma\chi_J)\,dH(y).                              \tag{2}
\]

Changing from counting measure to Haar probability only changes the
displayed kernel density convention, not (2), Markovness, or the estimates
below.

In global physical centroid coordinates the local cost is

\[
 (\mathcal J_cG)(u,X,\sigma)
 =\int_0^\infty e^{-s}K_c^s[G(u+s,\cdot)](X,\sigma)\,ds, \tag{3}
\]

on

\[
 d\nu=e^{-u}\,du\,dH(X)\,d\lambda(\sigma).              \tag{4}
\]

This is the complete finite-`c`, complete-label operator; no unlabelled or
unconditioned replacement is being made.

## Exact physical `L^p` norm

For `p>1`, the norm is exactly

\[
 \boxed{\ \|\mathcal J_c\|_{L^p(\nu)\to L^p(\nu)}={p\over p-1}.\ } \tag{5}
\]

For the upper bound, a fixed `s` summand in (3) has norm at most

\[
 e^{-s(1-1/p)}\|G\|_p,
\]

because `K_c^s` is angular-label Markov and `u\mapsto u+s` changes the
weight in (4) by `e^s`.  Minkowski and
`\int_0^\infty e^{-s/p'}ds=p'` give (5) with `<=`.

For completeness, the sharp lower bound does not require a limiting action.
Let `B_R` be the radius-`R` ball in `H_N`, take all labels, and define

\[
 G_{L,R}(u,X,\sigma)=e^{u/p}\mathbf1_{[0,L]}(u)
 \mathbf1_{B_R}(X).                                       \tag{6}
\]

If `S=s<=T`, then `\|PZ\|_2\le\|Z\|_2\le s`.  Hence for
`0<=u<=L-T` and `X\in B_{R-T}`, every summand with `0<=s<=T` in (3)
sees the value `e^{(u+s)/p}` in (6).  Consequently

\[
 {\|\mathcal J_cG_{L,R}\|_p\over\|G_{L,R}\|_p}
 \ge \left({L-T\over L}{|B_{R-T}|\over|B_R}\right)^{1/p}
 \int_0^T e^{-s/p'}ds.                                    \tag{7}
\]

First let `R,L\to\infty` for fixed `T`, then let `T\to\infty`.
This proves the reverse inequality in (5), uniformly in `c` and with the
full label system.  In particular the `p'` blow-up is an exact `L^p` fact,
not merely a loose Minkowski artifact.

## A sharp atom/column weak test, but not a weak-type theorem

The full finite-action kernel has a useful one-column calculation.  Normalize
a point source at `w_0=0` by its physical source mass.  With the Haar
probability convention fixed above, at radial separation `s` its field height
is

\[
 a_{\rm Haar}(s)={2^N e^{-s}\over\sqrt N\,H_{N,c}(s)},    \tag{8}
\]

and the aggregate receiver-shell measure, relative to the physical mass of
the source, is

\[
 e^s{\sqrt N\,H_{N,c}(s)\over2^N}\,ds.                    \tag{9}
\]

Indeed, (9) combines `e^{-u}/e^{-w}=e^s` for `w=u+s` with
`\sum_J|D_{s,J}|=\sqrt N H_{N,c}(s)` and the Haar mass `2^{-N}` of each
label.  The function `H_{N,c}` is increasing on `(0,\infty)` (strictly so
there), hence `a_{\rm Haar}` is decreasing.  For a level
`a_{\rm Haar}(q)`, (8)--(9) give

\[
 \begin{aligned}
 {\nu\{\hbox{column}>a_{\rm Haar}(q)\}\over\hbox{source mass}}
 &\le \int_0^q e^s{\sqrt N H_{N,c}(s)\over2^N}\,ds\\
 &\le {e^q\sqrt N H_{N,c}(q)\over2^N}
 ={1\over a_{\rm Haar}(q)}.                              \tag{10}
 \end{aligned}
\]

With counting measure on labels, the two factors `2^N` in (8)--(9) are both
absent, giving the same product.  Moreover this atom bound is sharp: since
`H_{N,c}` is eventually a positive polynomial of degree `N-1`, the change of
variables `t=q-s` gives

\[
 {\int_0^q e^sH_{N,c}(s)\,ds\over e^qH_{N,c}(q)}
 =\int_0^q e^{-t}{H_{N,c}(q-t)\over H_{N,c}(q)}\,dt\longrightarrow1.
\]

Dominated convergence applies because the quotient is at most one.  Thus
every individual column has weak quasi-norm exactly one.  This does
**not** extend by positivity to arbitrary superpositions of columns: that
is exactly the missing weak `(1,1)` issue.

## The `c=infty` physical kernel and its beta mixture

With no Bernoulli labels, `H_{N,\infty}(s)=s^{N-1}/(N-1)!`.  On one
physical orthant put `t_i=e^{-v_i}`.  The physical measure becomes ordinary
Lebesgue measure in `t` (up to the fixed coarea normalization), and (3) is
the explicit multiplicative operator

\[
 (\mathcal J_\infty F)(t)
 =\int_{(0,1)^N}{(N-1)!\over[-\log(r_1\cdots r_N)]^{N-1}}
 F(t_1r_1,\ldots,t_Nr_N)\,dr.                              \tag{11}
\]

For `\lambda>=1`, define the one-dimensional Hardy average

\[
 H_\lambda f(t)=\int_0^1\lambda r^{\lambda-1}f(tr)\,dr.
\]

Then the exact positive mixture is

\[
 \boxed{\quad
 \mathcal J_\infty=
 \int_1^\infty H_\lambda^{\otimes N}
 {(N-1)(\lambda-1)^{N-2}\over\lambda^N}\,d\lambda .\quad} \tag{12}
\]

Indeed, after multiplying the product kernels, the `\lambda` integral is
`(N-1)\int_0^\infty a^{N-2}e^{-as}da=(N-1)!/s^{N-1}` for
`s=-\log\prod r_i`.  Formula (12) explains both the exact `L^p` Hardy
constant and why a naive bound by the individual product Hardy weak norms is
unhelpful.

For `N=2`, the reciprocal Laplace/Fourier symbol can be written using the
logarithmic mean

\[
 L(a,b)={a-b\over\log a-\log b}=\int_0^1a^\theta b^{1-\theta}\,d\theta,
\]

so the symbol of `\mathcal J_\infty` is the associated logarithmic divided
difference (the reciprocal of `L` in this parametrization).  The integral
representation is suggestive, but no valid weak `(1,1)` deduction from it
has been supplied here.  In particular its positive mixed derivative does
not give the usual multivariable maximum-principle/Markov-generator route.

## Why no Hopf semigroup is presently available

For the unlabelled uniform-simplex conditional law, write `Z_s=sU` with
`U\sim\operatorname{Dirichlet}(1,\ldots,1)`.  Then

\[
 \operatorname{Cov}(PZ_s)={s^2\over N(N+1)}P.
\]

Independent conditional increments at `s,t` instead have covariance

\[
 {s^2+t^2\over N(N+1)}P,
\]

whereas the conditional increment at `s+t` has covariance
`(s+t)^2P/[N(N+1)]`.  Thus `K_sK_t\ne K_{s+t}` already at `c=infty`.
For finite c>0 choose s,t>0 with s+t<c. Every conditional Bernoulli label is then zero and exactly the same uniform-simplex covariance calculation applies. Thus the failure of this additive-parameter semigroup law also holds at the actual finite reference actions.

The usual Hopf/Dunford--Schwartz semigroup weak-endpoint theorem therefore
does not apply automatically.

## Consequence for the screening target

The original target `\|\mathcal J_c\|_{1\to1,\infty}=N^{o(1)}` is stronger
than the fixed-reference band screen needs.  If a weak bound with loss
`N^{\gamma+o(1)}` is proved for the actual fixed reference action, the
r069 screen has leading threshold `\gamma\le2\kappa`; at the endpoint take
`K=\kappa+\delta`, `\varepsilon=2\kappa+\delta`.  For the shifted reference,
`2\alpha=(9-\mu)/12\approx0.6434613`, so a hypothetical
`O(\sqrt N\log N)` weak bound would suffice.  The exact `L^p` result (5)
only yields an `O(N)` capped-screen coefficient and is not evidence for a
matching weak lower bound.  Establishing any exponent below `2\kappa` for
the full finite-action labelled operator remains open in this note.


The boundary exponent gamma=2*kappa also recovers the leading logarithmic coefficient by choosing K=kappa+delta and epsilon=2*kappa+delta for each fixed delta>0, then taking N to infinity before delta decreases to zero; screening70.md records this ordering.
