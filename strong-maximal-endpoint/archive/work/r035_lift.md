# r035 — tensor lifting the r023 near-unit corner gadget

Date: 2026-11-06.  This note examines only direct tensor/product lifts of
the r023 two-dimensional corner gadget.  It distinguishes an exact
Orlicz-denominator calculation from actual strong-rectangle product
witnesses.  It does not bound the full halo of a hyperbolic base input.

## The exact tensor denominator

Normalize the r023 gadget to the unit square:

\[
 g_\varepsilon=1+\varepsilon a,
 \qquad a=1_{[0,1]^2\setminus[0,1/4]^2},
 \qquad \alpha=\mathbb Ea=15/16.
\]

For `k` independent two-coordinate blocks, let
\(G_{k,\varepsilon}=\prod_{i=1}^k g_\varepsilon\).  On the unit
`2k`-cube, the number `J` of high gadget blocks has the exact
`Binomial(k,alpha)` law and

\[
 G_{k,\varepsilon}=(1+\varepsilon)^J.\tag{1}
\]

Consequently, if the remaining coordinates carry a constant amplitude
`A>=1`, the *actual* `Phi_n` denominator (not a product of `Phi_2`
denominators) is

\[
 \begin{aligned}
 D_{n,k}(A,\varepsilon)
 &=A\sum_{j=0}^k\binom{k}{j}\alpha^j(1-\alpha)^{k-j}(1+\varepsilon)^j\\
 &\quad\times[1+\log A+j\log(1+\varepsilon)]^{n-1}.
 \end{aligned}\tag{2}
\]

This formula applies whenever the input is `A G` on the full unit cube;
it includes the `n-2k` baseline cube coordinates, which contribute no
random factor.  Direct Wolfram 15.0 symbolic differentiation of (2) gives

\[
 \frac{\partial_\varepsilon D_{n,k}(A,0)}
  =A(1+\log A)^{n-1}\,k\alpha
    \left(1+\frac{n-1}{1+\log A}\right).\tag{3}
\]

For a bounded baseline, say `1<=A<=A_0`, the relative first-order
denominator price is therefore `Theta(n k epsilon)`, not the sum of the
two-dimensional prices `Theta(k epsilon)`.

The same fact holds for a separated hyperbolic block.  If `H(y)>=1` is a
bounded function of the other `n-2k` coordinates and
\(F(y,x)=H(y)G_{k,\varepsilon}(x)\), then differentiation under the
integral gives

\[
 \left.\frac d{d\varepsilon}\int\Phi_n(F)\right|_{0}
 =k\alpha\int\Phi_n(H)
 \left(1+\frac{n-1}{1+\log H}\right).\tag{4}
\]

In particular, for the cap-two profiles used in the fixed-block r026
experiments, the multiplier in parentheses is at least
\(1+(n-1)/(1+\log2)\), hence remains a fixed positive multiple of `n`.

## What a direct product witness can and cannot retain

For disjoint coordinate blocks, the actual strong operator factors exactly:

\[
 M_n\!\left(H\prod_{i=1}^k g_i\right)
 =M_{n-2k}H\prod_{i=1}^k M_2g_i.\tag{5}
\]

Thus one may take the Cartesian product of the r023 strict two-dimensional
witness regions (and a strict base witness).  Each gadget supplies only a
`1+O(epsilon)` relative-volume factor, so this particular simultaneous
product witness has volume multiplier `exp(O(k epsilon))`.  By (3)--(4),
its bounded-base Orlicz denominator has first-order multiplier
`exp(Theta(n k epsilon))`.  Therefore the r023 per-gadget first-order gain
does **not** survive this naive tensor/product-witness lift at high
dimension.

This is a scoped mechanism obstacle, not an upper bound on the full strong
halo: a receiver set that uses different subsets of active gadgets, or a
base profile tuned just below threshold, is not controlled by this argument.

## Full-cube special case: a rigorous ceiling

When there is no hyperbolic block and `A=1`, the tensor input is supported
on `Q=[0,1]^n` and satisfies

\[
 1\le G_{k,\varepsilon}\le(1+\varepsilon)^k\quad\hbox{on }Q.
\]

The r024 full near-flat-cube theorem therefore applies directly:

\[
 \mathcal L(G_{k,\varepsilon})
 \le V_n((1+\varepsilon)^k),\tag{6}
\]

where `V_n` is the exact constant-cube halo.  In particular, if
`k epsilon=o(1)`,

\[
 \log\mathcal L(G_{k,\varepsilon})
 \le O\!\big(\sqrt{n k\varepsilon}\big)+o(1),\tag{7}
\]

using the r024 estimate
\(V_n(1+t)\le(1+t)\exp(2\sqrt{2n\log(1+t)})\).  This rules out reading
the fixed-dimensional r023 first-order advantage as an additive `k gamma`
gain in the full-cube tensor family.  It still does not rule out thin
hyperbolic support, which lies outside the r024 class.

## Fixed-block cone comparison does not close the remaining case

The r026 conditional-density lemma for a single `m`-block corner gives the
factor

\[
 C_m\sqrt{\frac{n}{\prod_{j=1}^m b_j}}\,Q_n
\]

for a block partition of sizes `b_j`.  If the proposed gadgets are literally
`k` two-coordinate blocks, this becomes
\(C_k\sqrt{n/2^k}\,Q_n\).  It is a strong loss for balanced macroscopic
blocks, but it is *not* a uniform obstruction for fixed or slowly growing
two-dimensional blocks (and `C_k` was not controlled uniformly).  Thus it
would be incorrect to cite r026 alone as a no-go theorem for all tensor
lifts of r023.

## Result

The exact `Phi_n` coupling supplies a rigorous obstacle for direct bounded
base/product-witness lifts: the cost is order `n k epsilon`, whereas the
visible r023 witness-volume benefit is order `k epsilon`.  The pure-cube
tensor case additionally falls under the r024 ceiling.  A possible route
not excluded here would have to use thin support plus a non-product union of
receivers (or a large-log-amplitude base), and must calculate its real
`Phi_n` denominator rather than multiply two-dimensional ratios.

## Finite theorem for selected threshold-one Cartesian witnesses

There is a sharper, completely finite statement for the pure product input.
Let `n>=8`, `k<=n/2`, `0<epsilon<=1`, put
\(u=\log(1+\varepsilon)\) and \(x=ku\), and let
\(G=\prod_{i=1}^k g_\varepsilon\) (with unused coordinates equal to the
unit-cube baseline).  Select any Cartesian witness formed by taking, in
every active pair, points in the threshold-one halo
\(\{M_2g_\varepsilon>1\}\).  Then its volume divided by the **actual**
`Phi_n` denominator of `G` is at most one.

Indeed, since \(g_\varepsilon\le(1+\varepsilon)1_{[0,1]^2}\),

\[
 |\{M_2g_\varepsilon>1\}|
 \le V_2(1+\varepsilon)=1+4e^u u\le e^{4u}.\tag{8}
\]

The last inequality follows from differentiating
\(e^{4u}-1-4e^u u\).  Its derivative is
\(4e^u(e^{3u}-1-u)\ge0\).  Thus the selected Cartesian volume is at most
\(e^{4u})^k=e^{4x}\).

On the other hand, with `J~Binomial(k,15/16)`, the denominator is exactly

\[
 D=\mathbb E\big[e^{uJ}(1+uJ)^{n-1}\big].\tag{9}
\]

The integrand is convex in `J` (its second derivative is nonnegative for
`n>=2`), so Jensen gives

\[
 D\ge e^{\alpha x}(1+\alpha x)^{n-1},\qquad \alpha=15/16.\tag{10}
\]

It remains only to check

\[
 F_n(x):=(n-1)\log(1+\alpha x)-(4-\alpha)x\ge0,
 \quad 0\le x\le\frac n2\log2.\tag{11}
\]

`F_n` is concave in `x`, and `F_n(0)=0`, so it suffices to check its other
endpoint.  Write

\[
c=\frac{15}{32}\log2,\qquad d=\frac{49}{32}\log2,
\qquad F_n(n\log2/2)=(n-1)\log(1+cn)-dn.
\]

This endpoint expression is increasing for `n>=8`.  A fully rational
check uses

\[
\frac23<\log2<\frac7{10},\quad \frac5{16}<c<\frac{21}{64},\quad
d<\frac{343}{320},\quad \log(7/2)>\frac54.
\]

The last inequality follows from the first four positive terms of
\(2\operatorname{arctanh}(5/9)\).  At `n=8`, the endpoint is bounded
below by
\(7\cdot5/4-8\cdot343/320=7/40>0\); its derivative is bounded below by
\(1+35/58-343/320=4933/9280>0\).  This proves (11), and hence

\[
 \frac{|\text{selected Cartesian witness}|}{\int\Phi_n(G)}\le1.\tag{12}
\]

The scope of (12) is deliberately narrow: it bounds only the simultaneous
Cartesian product of the individual two-dimensional threshold-one halos.
It is not an upper bound for the full halo, does not cover receiver unions
with varying active gadget sets, and says nothing about a thin hyperbolic
base.

The exact Wolfram input and output for all rational inequalities are
`r035_naive_product_check.wl` and `r035_naive_product_check.json` in this
same external work directory.
