# Physical \(L^p\) screening and the weak-endpoint criterion

## The local-cost operator

Fix an action \(c\), retain the complete label group, and write
\[
 j_cG(u,X,\sigma)
 =\int_0^\infty e^{-s}
 [K_c^sG_{u+s}](X,\sigma)\,ds.
\]
Both source and receiver use the physical measure
\[
 d\mathfrak m=e^{-u}\,du\,dH(X)\,d\lambda(\sigma).
\]
Here \(K_c^s\) is a positive conditional joint-label Markov convolution,
so it is an \(L^p(dH\,d\lambda)\) contraction for every \(1\le p\le
\infty\).  A radial translation gives
\[
 \|e^{-s}K_c^sG_{u+s}\|_{L^p(d\mathfrak m)}
 \le e^{-s(1-1/p)}\|G\|_{L^p(d\mathfrak m)}.
\]
Minkowski therefore proves, for every \(p>1\),
\[
 \boxed{\ \|j_cG\|_p\le {p\over p-1}\|G\|_p.\ }             \tag{1}
\]
This is uniform in \(c\), including the limiting action, and includes
all labels.

The constant \(p'=p/(p-1)\) is sharp for this argument and for the
operator norm: use a large angular box on which conditional translations
have negligible boundary loss, and the usual long, truncated radial Hardy
near-extremizer.  Thus no improvement of (1) to a bounded \(p\downarrow1\)
constant can be expected.

## What the \(L^p\) bound actually screens

Let \(0\le G\le M\), put \(J(G)=\int G\,d\mathfrak m\), and take
\[
 H=\log(M/t)\ge1,\qquad p=1+H^{-1}.
\]
Since \(\|G\|_p^p\le M^{p-1}J(G)\), Chebyshev and (1) give
\[
 \begin{aligned}
 \mathfrak m\{j_cG>t\}
 &\le t^{-p}(p')^pM^{p-1}J(G)\\
 &={J(G)\over t}\,e\,(H+1)^{1+1/H}\\
 &\le {2e(H+1)J(G)\over t}.                                \tag{2}
 \end{aligned}
\]
This is an \(L\log L\)-type screening estimate, not weak type \((1,1)\).
No different single choice of \(p\) repairs its polynomial loss.  The exact
best constant obtainable from (1) and Chebyshev is
\[
 C(H)=\inf_{p>1}(p')^p e^{(p-1)H}.
\]
Writing \(q=p'>1\), its unique stationary point satisfies
\[
 q-\log q=H+1.
\]
It is the minimum, and
\[
 q=H+\log H+1+o(\log H),\qquad C(H)=eq\sim eH.
\]
Thus the \(\Theta(N)\) loss below is a one-shot
\(L^p\)-Chebyshev method obstruction, not a lower bound on the actual weak
norm.

For \(G=\Phi_N(F)\), \(F\le A\) with fixed \(A>1\), one has
\(\log M=\Theta_A(N)\).  If
\[
 t=J_{\rm ref}N^\varepsilon,
\]
then \(H=\Theta_A(N)\), and (2) only gives
\[
 \mathfrak m\{j_cG>t\}
 =O_A\!\left(Q_{\rm ref}N^{1-\varepsilon}J(G)\right).       \tag{3}
\]
In particular, for \(J(G)\le J_{\rm ref}=Q_{\rm ref}^{-1}\), its useful
decay requires \(\varepsilon>1\).  Formula (3) is only a bound on the
discarded screening set.  It is neither a lower bound on an actual halo
nor a proof of a weak-\((1,1)\) estimate.  The boundedness
\(G\le\Phi_N(A)\) is essential here; this screen does not cover uncapped
sources.

## Compatibility with forbidden-band capture

The scalar forbidden-band result in r069 says that a source admitted at
budget \(J_{\rm ref}N^\varepsilon\) can be forced into a band with
parameter \(K\) only under
\[
 \varepsilon<2K,
\]
and its retained band budget is lower by
\[
 N^{K-\kappa+o(1)}
\]
relative to the calibrated reference scale.  For the near-optimal choice
\(K=\kappa+\delta\), the \(L^p\) screen would require both
\(\varepsilon>1\) and \(\varepsilon<2\kappa+o(1)\).  This is impossible
for the shifted scalar reference because
\[
 2\alpha={9-\mu\over12}<1.
\]
The fixed-\(c\) coefficient is smaller still.  Thus (1) cannot supply the
near-logarithmic-scale screening needed by that bridge.

## A sufficient weak-endpoint exponent

Suppose instead that one proves, for a fixed reference action and all
nonnegative \(G\), the physical weak estimate
\[
 \mathfrak m\{j_cG>t\}
 \le C N^{\gamma+o(1)}{J(G)\over t}.                        \tag{4}
\]
At \(t=J_{\rm ref}N^\varepsilon\), this makes the bad screening set
\(N^{\gamma-\varepsilon+o(1)}\) times the normalized source budget.
It is compatible with r069 whenever one can choose
\[
 K>\max\{\kappa,\gamma/2\},\qquad \gamma<\varepsilon<2K.    \tag{5}
\]
Letting the strict gaps in (5) tend to zero gives the qualified
fixed-reference endpoint exponent
\[
 \kappa+\max\{0,\gamma/2-\kappa\}
 =\max\{\kappa,\gamma/2\}.                                  \tag{6}
\]
Equivalently, a weak loss \(\gamma<2\kappa\) is already sufficient for
the near-\(\kappa\) scalar scale; a proof of \(N^{o(1)}\) weak type is not
necessary.  The boundary \(\gamma=2\kappa\) has the same leading
coefficient: for every \(\delta>0\), choose
\[
 K=\kappa+\delta,\qquad \varepsilon=2\kappa+\delta.
\]
Then \(\varepsilon<2K\), the bad-set exponent is \(-\delta\), and the
capture loss is only \(\delta\); take \(\delta\downarrow0\) after the
asymptotic statement.  Thus \(\gamma<2\kappa\) is a convenient strict
sufficient condition, while the leading-coefficient threshold is
\(\gamma\le2\kappa\), in agreement with (6).

The transfer from (4)--(6) to an actual reference halo still requires the
separate local ingredients: reference mean \(1-o(L/N)\), a common-source
comparison, and a packing/stopping step.  Nothing here establishes those
ingredients or a general variable-action upper bound.  Fixed \(c\)
real-average results such as r028 are stronger in their own setting; this
note records only the bridge criterion.
