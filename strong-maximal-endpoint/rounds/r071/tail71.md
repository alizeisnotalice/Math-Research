# Fixed-dimensional tails, critical columns, and posterior truncation

## Scope

Assume throughout that G is nonnegative. All assertions in the first section are for fixed \(N\ge2\), fixed
\(c>0\), and a source supported in the logarithmic Cartesian box
\(y_i\in[-R,R]\).  They are not uniform in the critical regime
\(N\to\infty\).

Use full counted labels and write
\[
 m=\sum_\sigma\int e^{-\sum_i y_i}G_\sigma(y)\,dy<\infty.
\]
Let \(j_cG\) be the full local-cost field in log-Cartesian coordinates.
If \(m=0\), then \(G=0\) almost everywhere and all tail assertions below
are immediate.  Hence assume \(m>0\).

## Exact fixed-\(N\) weak tail

Here \(H_{N,c}(s)=\Gamma(N)^{-1}\sum_{k=0}^N{N\choose k}(s-ck)_+^{N-1}\). The Cartesian kernel for a fixed receiver label is \(e^{-S}/H_{N,c}(S)\), integrated over source increments \(z_i\ge cJ_i\), summed over increment labels \(J\), with source label \(\sigma\mathbin\oplus J\). Thus the field in (2) is bounded by \(e^{-T}\sum_\tau\int e^{-\sum y}G_\tau(y)/H_{N,c}(T+\sum y)dy\).

Put \(T=-\sum_i x_i\).  The joint-label normalizer has, uniformly for
\(|w|\le NR\),
\[
 H_{N,c}(T+w)
 ={2^N\over\Gamma(N)}T^{N-1}
 \left(1+O_{N,c,R}(T^{-1})\right).                          \tag{1}
\]
The cancellation of the exponential source density against the physical
Jacobian and the complete-label conditional density therefore gives the
uniform upper field estimate
\[
 j_cG(x,\sigma)\le a(T)(1+O_{N,c,R}(T^{-1})),\qquad
 a(T)={m\Gamma(N)e^{-T}\over2^NT^{N-1}},                    \tag{2}
\]
whenever a receiver can be reached from the source.  Such a receiver
satisfies \(x_i\le R\).  On the outer cone
\[
 x_i\le-R-c\quad(1\le i\le N),
\]
every source point and every required label is admissible, and the same
calculation gives the matching lower estimate in (2).

The counted-label cumulative physical cone volumes at level T, between the corresponding
upper and lower cones are
\[
 {2^Ne^T\over\Gamma(N)}
 (T\pm N(R+c))^{N-1}(1+O_{N,R}(T^{-1})).                    \tag{3}
\]
For a small level \(t\), let \(q=q(t)\to\infty\) solve
\[
 a(q)=t.
\]
Since
\[
 -\partial_T\log a(T)=1+{N-1\over T}\ge1,
\]
the relative \(O(T^{-1})\) field errors in (2) move the level cutoff only
by \(O_{N,c,R}(q^{-1})\), not by \(O(1)\).  More explicitly, for a fixed
constant \(C\) and all large \(q\), the lower outer cone cut at
\(T=q-C/q\) is contained in \(\{j_cG>t\}\), while the latter is contained
in the reachable upper cone cut at \(T=q+C/q\).  The cumulative cone
volumes are given by (3).  Multiplication by
\[
 t={m\Gamma(N)e^{-q}\over2^Nq^{N-1}}
\]
gives on both sides
\[
 m\,e^{\pm C/q}
 \left(1\pm{C\over q^2}\right)^{N-1}(1+O_{N,c,R}(q^{-1}))
 =m(1+O_{N,c,R}(q^{-1})).
\]
The squeeze therefore gives
\[
 \boxed{\displaystyle
 \lim_{t\downarrow0}t\,\bigl|\{j_cG>t\}\bigr|_{\rm physical}=m.} \tag{4}
\]
The receiver set with \(T\) bounded above has finite physical measure:
reachability gives \(x_i\le R\), while
\(T\le T_0\) gives \(x_i\ge-T_0-(N-1)R\).  It does not affect (4).

This proves a fixed-\(N\) tail identity, not a dimension-uniform
weak-\((1,1)\) estimate.

## A noncommuting critical-column limit

For one critical column define
\[
 D_N(s)=\int_0^s e^{-r}{H_{N,c}(s-r)\over H_{N,c}(s)}\,dr.  \tag{5}
\]
For fixed \(N,c\), \(H_{N,c}\) is a polynomial eventually, so for every
fixed \(r\),
\[
 {H_{N,c}(s-r)\over H_{N,c}(s)}\longrightarrow1,
\qquad D_N(s)\longrightarrow1\quad(s\to\infty).             \tag{6}
\]

In contrast, in the critical window
\(c=\mu+O(\sqrt{\log N/N})\), \(s=q_\circ+O(\log N)\), the uniform rate
estimate gives
\[
 \partial_s\log H_{N,c}(s)=1+o(1).
\]
For fixed \(r\), the ratio in (5) tends to \(e^{-r}\).  Since \(H\) is
increasing, it is at most one, and dominated convergence gives
\[
 D_N(s)\longrightarrow\int_0^\infty e^{-2r}\,dr={1\over2}. \tag{7}
\]
Thus the fixed-\(N\), deep-tail value one cannot be substituted for the
critical-window value one half.  The two limits do not commute.

### Quantitative compact-column refinement

Suppose \(c\) and \(s/N\) range in fixed positive compact sets, and that
\([s-\sqrt N,s]\) remains in the same macro-compact regime.  The uniform
rate formula from r068 is
\[
 B_{N,c}(s):=\partial_s\log H_{N,c}(s)
 =\gamma(s)+O(N^{-1}),
\]
where the tilted mean equation gives
\[
 \gamma'(s)=-{1\over N\,\operatorname {Var}_{\rm tilt}(Z)}
 =O(N^{-1}).                                                \tag{8}
\]
The \(O(N^{-1})\) remainder need not be differentiated: integrating its
two endpoint values contributes \(O(t/N)\), while (8) contributes
\(O(t^2/N)\).  Uniformly for \(0\le t\le\sqrt N\),
\[
 \int_0^tB_{N,c}(s-v)\,dv
 =B_{N,c}(s)t+O((t+t^2)/N).                                 \tag{9}
\]
Consequently, on this range,
\[
 {H_{N,c}(s-t)\over H_{N,c}(s)}
 =\exp\{-B_{N,c}(s)t+O((t+t^2)/N)\}.
\]
The part of (5) with \(t>\sqrt N\) is at most \(e^{-\sqrt N}\), since
\(H\) is increasing.  Integrating the preceding error against the
exponentially decreasing \(e^{-t}\) proves
\[
 \boxed{\displaystyle
 D_N(s)={1\over1+B_{N,c}(s)}+O(N^{-1})
 ={1\over1+\gamma(s)}+O(N^{-1}).}                           \tag{10}
\]

At \(c=\mu\), \(s=q_\circ\), the mean equation expanded at the untilted
point gives
\[
 \gamma(q_\circ)
 =1+{z_N\over\sqrt{\mu N}}+O(\log N/N).
\]
Hence
\[
 D_N(q_\circ)
 ={1\over2}-{z_N\over4\sqrt{\mu N}}+O(\log N/N).            \tag{11}
\]
This is a conditional-column calculation.  Its column level
\(a_N(s)\) is not the r069 screening threshold
\(J_{\rm ref}N^\varepsilon\); (10)--(11) therefore do not by themselves
give a new bound for \(C_N\).

## A posterior truncation lemma

For fixed \(s>0,c>0\), let
\[
 w_k={N\choose k}(s-ck)_+^{N-1},\qquad k\ge0,
\]
with ratios on their positive support
\[
 r_k={w_{k+1}\over w_k}
 ={N-k\over k+1}
 \left({s-c(k+1)\over s-ck}\right)^{N-1}.                  \tag{12}
\]
The concavity of \(k\mapsto\log(s-ck)\) and the binomial factor give
\[
 \log r_k-\log r_{k-1}
 \le\log {k(N-k)\over(k+1)(N-k+1)}
 \le-{4\over N+2}.                                         \tag{13}
\]
For the final inequality, write the logarithm as
\[
 \log\left(1-{1\over k+1}\right)
 +\log\left(1-{1\over N-k+1}\right)
\]
and use \(\log(1-z)\le-z\) together with
\((k+1)^{-1}+(N-k+1)^{-1}\ge4/(N+2)\).

Choose any mode \(m\), so that \(r_m\le1\le r_{m-1}\), with the evident
one-sided interpretation at an endpoint.  Iterating (13) gives
\[
 {w_{m\pm j}\over w_m}
 \le\exp\left\{-{2j(j-1)\over N+2}\right\},\qquad j\ge1,    \tag{14}
\]
where nonexistent indices are omitted.  Let
\[
 R_{\rm cut}=\left\lceil\sqrt{110(N+2)}\right\rceil.
\]
Every omitted index from the truncation \(|k-m|\le R_{\rm cut}\) has
\(j\ge R_{\rm cut}+1\), hence its bound in (14) is at most \(e^{-220}\).
Since \(\sum_k w_k\ge w_m\),
\[
 {\sum_{|k-m|>R_{\rm cut}}w_k\over\sum_kw_k}
 \le (N+1)e^{-220}.                                       \tag{15}
\]
Consequently replacing the exact posterior by the normalized retained
posterior changes the expectation of every test in \([0,1]\) by at most
\((N+1)e^{-220}\).  This is an exact finite-\(N\) deterministic
truncation bound.

The posterior is over the count k. If used for a labelled observable, retain its exact conditional label distribution at each k; truncating k does not authorize arbitrary label averaging.
