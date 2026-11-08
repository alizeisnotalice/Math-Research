# r081: soft boundary density for fixed-corner source columns

## 1. A compact \(d=N-1\) dimensional boundary chart

Let
\[
 F_N(s)=2\prod_{i=1}^Ns_i\,\mathbf1_{[0,1]^N}(s),\qquad
 u(x)=x-\sqrt{x^2-1},\qquad
 X(v)={v+v^{-1}\over2}.
\]
For \(x_i>1\), the one-dimensional maximization from r080 tensorizes:
\[
 M_{\rm out}F_N(x)=2\prod_{i=1}^Nu(x_i),                 \tag{1}
\]
and the unique maximizing corner is \(a_i=u(x_i)\).

Put \(d=N-1\) and \(q=2^{-1/N}\).  Choose
\[
 0<a<q<b<1,
\]
then a compact box \(K\Subset(a,b)^d\) containing
\((q,\ldots,q)\), and finally \(\delta_0>0\), so small that
\[
 v_N(v,\delta):={1+\delta\over2\prod_{i=1}^dv_i}\in[a,b]
 \quad(v\in K,\ 0\leq\delta\leq\delta_0).               \tag{2}
\]
The compact receiver patch is parameterized by
\[
 T(v,\delta)=\bigl(X(v_1),\ldots,X(v_d),
                    X(v_N(v,\delta))\bigr).             \tag{3}
\]
It has depth \(M_{\rm out}F_N(T(v,\delta))-1=\delta\).
The map is a smooth diffeomorphism onto its image.  Indeed, on \([a,b]\),
\[
 \ell:={b^{-2}-1\over2}\leq|X'|\leq
 L:={a^{-2}-1\over2},
\quad {1\over2b^d}\leq\partial_\delta v_N\leq{1\over2a^d},
\]
so
\[
 {\ell^N\over2b^d}\leq|\det DT|\leq {L^N\over2a^d}.     \tag{4}
\]
All constants below may depend on this fixed compact chart and on \(N\).

## 2. One fixed true column has a square-root section

For \(v\in[a,b]\) and \(c\in[0,1]\), write
\[
 r_c(v)={1-c^2\over1+v^2-2cv}
 =1-{(c-v)^2\over1+v^2-2cv}\leq1.                       \tag{5}
\]
At (3), the anchored average of \(F_N\) equals
\[
 (1+\delta)\prod_{i=1}^N r_{c_i}(v_i).                  \tag{6}
\]
For sufficiently small \(\delta_0\), a successful column has every corner
coordinate in \((0,1)\).  A coordinate at least one sees no source in that
direction.  If one is at most zero, its one-dimensional factor is at most
\(1/(2x)\), whose ratio to \(u(x)\) is
\((1+v^2)^{-1}\leq(1+a^2)^{-1}\); taking \(\delta_0<a^2\) precludes
success.

If \(0\leq G\leq F_N\) succeeds from fixed corner \(c\) at (3), every
factor in (6) must exceed \((1+\delta)^{-1}\).  Since
\(1+v_i^2-2c_iv_i\leq1+b^2\),
\[
 |v_i-c_i|<\sqrt{1+b^2}\,\sqrt\delta.                   \tag{7}
\]
Hence its depth-\(\delta\) success section, in the \(v\in K\) variables,
has measure at most
\[
 C_{\rm sec}\delta^{d/2},\qquad
 C_{\rm sec}=\bigl(2\sqrt{1+b^2}\bigr)^d.                \tag{8}
\]

Put \(\alpha=X(b)-1>0\).  Every successful rectangle in this patch has
side lengths at least \(\alpha\), so
\[
 J_N(G)\geq\int G>|R|\geq\alpha^N=:c_N.                 \tag{9}
\]
For a countable weighted family set
\[
 \Theta(v,\delta)=\sum_jw_j\mathbf1_{\{G_j\ {\rm succeeds\ at}\ T(v,\delta)\}},
 \qquad B=\sum_jw_jJ_N(G_j).
\]
Then
\[
 \int_K\Theta(v,\delta)\,dv
 \leq {C_{\rm sec}\over c_N}B\,\delta^{d/2}
 =:C_N^\flat B\,\delta^{d/2}.                           \tag{10}
\]
This is only a fixed-corner source-column estimate.

## 3. Polynomial soft densities

Fix \(\lambda>0\).  A soft target of order \(\beta\geq0\) means
\[
 \Theta(v,\delta)\geq\lambda\delta^\beta
\quad\hbox{for a.e. }(v,\delta)\in K\times(0,\delta_0). \tag{11}
\]

### \(\beta<d/2\): finite budget fails, with positive deficit

For almost every depth, (10)--(11) imply
\[
 \lambda|K|\delta^\beta\leq C_N^\flat B\delta^{d/2}.     \tag{12}
\]
Thus finite \(B\) cannot give exact soft coverage.  Let
\[
 q_0=d/2-\beta>0,\qquad
 \delta_B=\min\left\{\delta_0,
 \left({\lambda|K|\over2C_N^\flat B}\right)^{1/q_0}\right\}.
\]
The soft deficit
\[
 D_\beta=\int_0^{\delta_0}\!\!\int_K
       (\lambda\delta^\beta-\Theta(v,\delta))_+\,dv\,d\delta
\]
satisfies
\[
 D_\beta\geq{\lambda|K|\over2(\beta+1)}\,
             \delta_B^{\beta+1}>0.                     \tag{13}
\]
Therefore a uniform finite budget cannot make this \(L^1\) soft deficit
vanish.

### \(\beta=d/2\): exact failure, but uniform \(L^1\) approximation

The slice integral is now scale-neutral.  Exact a.e. coverage is still
impossible for a countable finite-budget family.  Associate to it the
finite atomic measure
\[
 \mu=\sum_jw_j\delta_{(c_{j1},\ldots,c_{jd})}
\quad\hbox{on }\mathbb R^d.
\]
After discarding irrelevant columns, (9) gives
\(\mu(\mathbb R^d)\leq B/c_N\).  By (7), (11) would force, for almost every
\((v,\delta)\),
\[
 \mu(B_\infty(v,C\sqrt\delta))\geq\lambda\delta^{d/2}.
\tag{14}
\]

Here is an a.e.-valid tail contradiction.  Choose finitely many atoms so
that the omitted tail \(\nu\) has mass below a sufficiently small constant
multiple of \(\lambda|K|\).  For every sufficiently small \(t\), the union
of the \(Ct\)-balls around these finitely many atoms has measure less than
\(|K|/2\).  On the remaining set \(A_t\subset K\), Fubini gives
\[
 \int_{A_t}\nu(B_\infty(v,Ct))\,dv\leq C'_N t^d\nu(\mathbb R^d).
\]
Markov's inequality and the chosen tail mass show that a subset of
\(A_t\) of positive measure has
\(\mu(B_\infty(v,Ct))<\lambda t^d\).  Since (11) holds for almost every
pair, choose a sufficiently small \(\delta\) from its full-measure set and
put \(t=\sqrt\delta\).  This yields a positive-measure set of violations
of (14), a contradiction.  Thus finite atomic budget cannot give exact
critical soft coverage.

Uniform \(L^1\) approximation is different.  Extend \(K\) by a fixed
compact collar inside the open chart domain where (2) remains valid, and place an \(h\)-net there with weight
\(C\lambda h^d\) on each column.  Its total weight is \(O_N(\lambda)\).
For \(\delta\gtrsim h^2\), every \(v\in K\) has
\(\gtrsim(\sqrt\delta/h)^d\) net points within a \(c\sqrt\delta\)-ball,
including at the boundary of \(K\), so their combined weight is
\(\gtrsim\lambda\delta^{d/2}\).  The sufficient success estimate below
shows those columns really succeed.  The omitted layer has soft mass
\[
 O_N\!\left(\int_0^{Ch^2}\delta^{d/2}\,d\delta\right)
 =O_N(h^{d+2})\longrightarrow0.                         \tag{15}
\]

### \(\beta>d/2\): a finite-budget exact soft cover

The needed sufficient estimate is explicit.  On the compact chart,
\[
 1-r_c(v)\leq{|c-v|^2\over1-b^2}.                        \tag{16}
\]
Also \(v_N(v,\delta)\) is Lipschitz in \((v,\delta)\).  If a boundary
net corner is \(c^\ast=(v^\ast,v_N(v^\ast,0))\), then
\[
 |v-v^\ast|_\infty\leq c\sqrt\delta
\quad\Longrightarrow\quad
 |v_N(v,\delta)-c_N^\ast|\leq C_\ast(c\sqrt\delta+\delta).
\]
Choose \(c\) and then \(\delta_0\) so that the sum of the \(N\) losses in
(16) is at most \(\delta/2\).  Then
\[
 (1+\delta)\prod_{i=1}^Nr_{c_i^\ast}(v_i)
 \geq(1+\delta)(1-\delta/2)>1.                           \tag{17}
\]

Let \(\delta_k=2^{-k}\delta_0\).  For each \(k\), take a net of mesh
\(c\sqrt{\delta_{k+1}}\) in the collar and use the full source \(F_N\) at
each boundary-optimal net corner, with weight \(C\lambda\delta_k^\beta\).
By (17), the \(k\)-net covers the band
\(\delta_{k+1}<\delta\leq\delta_k\) with weight at least
\(\lambda\delta^\beta\).  Its total modular cost is bounded by
\[
 C_{N,\lambda}J_N(F_N)
 \sum_{k\geq0}\delta_k^{-d/2}\delta_k^\beta
 <\infty.                                                \tag{18}
\]

## 4. Flat density removes the fixed-\(N\) entropy obstruction

Let
\[
 \chi(z)={e^{-1/z}\over e^{-1/z}+e^{-1/(1-z)}}\quad(0<z<1).
\]
For any fixed \(\eta>0\), first shrink \(\delta_0\leq\eta\).  The same
dyadic construction with weights
\(C\lambda\chi(\delta_k/\eta)\) covers the target
\[
 \Theta(v,\delta)\geq\lambda\chi(\delta/\eta)
\]
and costs at most
\[
 C_{N,\lambda}J_N(F_N)
 \sum_{k\geq0}\delta_k^{-d/2}\chi(\delta_k/\eta)<\infty. \tag{19}
\]
This applies in particular to the current scale \(\eta=N^{-2}\).  The
constant in (19) may be extremely bad in \(N\) and \(\eta\); this removes a
polynomial entropy obstruction at each fixed dimension, not a
dimension-sharp budget barrier.

Nothing here implies a smooth-norm blowup or a strong-maximal constant
blowup.  It concerns only one global fixed-corner source-column budget
architecture.
