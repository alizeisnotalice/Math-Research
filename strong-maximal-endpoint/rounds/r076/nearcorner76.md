# Near-corner cap-two Hardy charge

## Statement

Let \(0\le F\le2\) be compactly supported.  Suppose each actual successful
rectangle \(R_b(x)=\prod_i[b_i(x),x_i]\), of one fixed orientation, is
contained in one member \(R_a(x)=\prod_i[a_i,x_i]\) of a prescribed finite
family of oriented fixed-corner rectangles, with
\[
 {|R_a(x)|\over|R_b(x)|}\le e^\tau.
\tag{1}
\]
Then actual success implies \(H_aF(x)>e^{-\tau}\).  The best cap-two
tensor-Hardy charge obtainable by this \(L^p\) argument is
\[
 \log B_N(\tau)=
 \inf_{1<p\le p_N}\left\{Np\log{p\over p-1}+p\tau\right\},
 \qquad p_N=1+\zeta(N-1),\quad
 \zeta={\log(1+\log2)\over\log2}.
\tag{2}
\]
For \(\tau=o(N)\), uniformly including \(\tau\asymp N^{-1}\),
\[
 \log B_N(\tau)=N+\sqrt{2N\tau}+O(\tau+1).
\tag{3}
\]

This is a finite-family conditional charge.  It does not construct a
continuous corner net for arbitrary policies.

## 1. Geometry and orientations

For nonnegative \(F\), (1) and containment give
\[
 H_aF(x)={1\over|R_a(x)|}\int_{R_a(x)}F
 \ge {|R_b(x)|\over|R_a(x)|}H_bF(x)>e^{-\tau}.
\tag{4}
\]
For one fixed anchor \(a\), the full signed-corner operator handles all
\(2^N\) receiver orientations without a \(2^N\) entropy loss.  After
translation and reflections, apply the one-sided tensor Hardy estimate on
each source orthant.  Those orthant source pieces are disjoint, so their
Orlicz modulars sum to at most the global modular.  Thus the same-anchor
orientation sum is already included in the one-anchor bound below.

An extra count is needed only for corners that cannot share this fixed
anchor/source-orthant partition (for example, genuinely different anchor
locations in a proposed cover).  The symbol \(M_N\) below denotes that
post-sharing count.

## 2. Cap-two Orlicz comparison and the exact minimizer

As in the fixed-corner calculation, for every \(1<p\le p_N\),
\[
 t^p\le\Phi_N(t)\qquad(0\le t\le2).
\tag{5}
\]
The tensor Hardy norm is \((p')^N\) on \(L^p\).  Chebyshev at height
\(e^{-\tau}\) therefore yields for one oriented corner
\[
 |\{H_aF>e^{-\tau}\}|
 \le (p')^{Np}e^{p\tau}\int\Phi_N(F).
\tag{6}
\]
Optimizing (6) gives (2).

Put \(q=(p-1)^{-1}\); the constraint \(p\le p_N\) becomes
\[
 q\ge q_{\min}:={1\over\zeta(N-1)}.
\tag{7}
\]
The objective in (2) is
\[
 L_N(q)=N(1+q^{-1})\log(1+q)+(1+q^{-1})\tau,
\tag{8}
\]
and direct differentiation gives
\[
 L_N'(q)={N[q-\log(1+q)]-\tau\over q^2}.
\tag{9}
\]
Thus the optimizer is:

\[
 \begin{cases}
 q=q_{\min},&0\le\tau\le\tau_{\rm crit},\\
 q>q_{\min}\text{ solving }q-\log(1+q)=\tau/N,
     &\tau>\tau_{\rm crit},
 \end{cases}
\qquad
 \tau_{\rm crit}=N[q_{\min}-\log(1+q_{\min})].
\tag{10}
\]
At an interior optimizer, substituting
\(\tau=N[q-\log(1+q)]\) into (8) gives the exact identity
\[
 \log B_N(\tau)=N(1+q).
\tag{11}
\]

## 3. Uniform small- and moderate-\(\tau\) asymptotics

Since \(q_{\min}=O(N^{-1})\), one has
\(\tau_{\rm crit}=O(N^{-1})\).  In the boundary regime of (10),
\[
 L_N(q_{\min})=N+{1\over2\zeta}+p_N\tau+O(N^{-1})=N+O(1),
\tag{12}
\]
which is consistent with (3), because \(\sqrt{N\tau}=O(1)\) there.

For the interior regime with \(u=\tau/N=o(1)\),
\[
 q-\log(1+q)={q^2\over2}+O(q^3),
 \qquad q=\sqrt{2u}+O(u).
\tag{13}
\]
Combining (11) and (13) gives
\[
 N(1+q)=N+\sqrt{2N\tau}+O(\tau),
\tag{14}
\]
which together with (12) proves the uniform remainder in (3), including
the transition \(\tau\asymp N^{-1}\).

## 4. Finite-corner entropy and the linear frontier

If the total number of post-sharing fixed-corner families is \(M_N\), let
\(m_N=\log M_N\).  The union bound and (6) give the charge
\[
 |\{\text{actual successes covered by the family}\}|
 \le M_NB_N(\tau)\int\Phi_N(F).
\tag{15}
\]
At the critical reference scale
\(\log Q_{\rm ref}=q_0+\alpha\log N+o(\log N)\) with
\(q_0/N\to\mu\), the leading criterion in the interior regime is
\[
 {m_N\over N}+1+q<\mu.                                  \tag{16}
\]
This compares the actual finite-family entropy and Hardy charge directly
with the endpoint budget; lower-order \(\sqrt{N\log N}\) and logarithmic
terms must still be retained at equality.

For \(\tau=\beta N\), the interior equation is
\[
 q-\log(1+q)=\beta.
\tag{17}
\]
With zero leading corner entropy, the frontier in (16) occurs at
\(q=\mu-1\), hence
\[
 \boxed{\ \beta_*=\mu-1-\log\mu\ }.                    \tag{18}
\]
For a genuinely finite or subexponential family, near-corner enlargements
below this linear threshold are globally chargeable by this method.

The full orientation sum at one fixed anchor is free by the disjoint-orthant
argument above.  This proof still does not claim that an arbitrary continuum
of moving anchors is automatically controlled: different anchors can require
the genuine entropy \(m_N\) in (16).  It only supplies the stated finite
post-sharing corner-family lemma for arbitrary cap-two sources.
