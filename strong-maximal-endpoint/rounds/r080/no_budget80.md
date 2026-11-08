# A two-dimensional obstruction to finite full-halo fixed-corner budgets

## Statement and scope

Let
\[
 F(s,t)=2st\,{\bf1}_{[0,1]^2}(s,t),
 \qquad \Phi_2(r)=r(1+\log^+r).                         \tag{1}
\]
For a real fixed corner \((a,b)\), write \(H_{a,b}G(x,y)\) for the
average of \(G\) over the oriented rectangle having \((a,b)\) and
\((x,y)\) as opposite corners, when it is nondegenerate.  A source column
is \((a,b,G,w)\), where \(0\le G\le F\), \(w\ge0\), and its success set is
\(\{H_{a,b}G>1\}\).  The charged cost is \(wJ(G)\),
\(J(G)=\int\Phi_2(G)\).

There is a finite-measure strict halo subset \(E\) such that no countable
fixed-corner source-column cover of \(E\) can have finite total cost.  More
quantitatively, every weighted family of columns of total cost \(B\) has
an \(L^1(E)\) coverage deficit bounded below by a positive function of
\(B\).  This rules out a proposed *full-halo fixed-corner source-column
budget* already in dimension two.  It does not contradict same-anchor
consolidation, and it is not a counterexample to a strong-maximal inequality
or a statement about \(C_N\).

## 1. The strict halo patch

For \(x>1\), set
\[
 u(x)=x-\sqrt{x^2-1}\in(0,1),\qquad
 K_a(x)={1-a^2\over2(x-a)}\quad(0\le a<1).             \tag{2}
\]
The identity
\[
 K_a(x)=u(x)-{(a-u(x))^2\over2(x-a)}                    \tag{3}
\]
shows that \(K_a(x)\le u(x)\), with equality uniquely at \(a=u(x)\).
For \(x,y>1\), taking the product of the one-dimensional maximizations
gives
\[
 M_{\rm S}F(x,y)=2u(x)u(y).                              \tag{4}
\]

Let
\[
 I=[25/24,13/12],\qquad |I|=1/24,
\tag{5}
\]
so \(u(I)=[2/3,3/4]\).  Define \(Y:I\to I\) by
\[
 u(Y(x))={1\over2u(x)}.                                  \tag{6}
\]
For \(0<\delta<\delta_0:=1/96\), put
\[
 \Gamma_\delta=\{(x,Y(x)-\delta):x\in I\},
 \qquad
 E=\{(x,Y(x)-\delta):x\in I,\ 0<\delta<\delta_0\}.   \tag{7}
\]
Since \(u\) decreases, (4) and (6) give
\[
 2u(x)u(Y(x)-\delta)>1.                                  \tag{8}
\]
Thus \(E\) is a strict subset of the genuine strong-maximal halo.  The
map \((x,\delta)\mapsto(x,Y(x)-\delta)\) has Jacobian one, so
\(|E|=|I|\delta_0<\infty\).

## 2. Every successful column has a fixed cost

A column successful at a point of \(E\) must have \(a,b\in(0,1)\).  Indeed,
an anchor coordinate at least one sees no positive-measure part of the
source in that direction.  If \(a\le0\), the first one-dimensional factor
is at most \(1/(2x)\), while the optimal second factor is at most \(u(y)\).
Hence
\[
 H_{a,b}F(x,y)\le {u(y)\over x}<1.                      \tag{9}
\]
The same argument applies if \(b\le0\).  Since \(G\le F\), success for
\(G\) implies success for \(F\), proving the assertion.

On \(E\), \(x-1\ge1/24\) and
\(y-1\ge33/32-1=1/32\).  For \(a,b\in(0,1)\), the successful rectangle
has area
\[
 (x-a)(y-b)\ge(x-1)(y-1)\ge1/768.                       \tag{10}
\]
If \(H_{a,b}G(x,y)>1\), then
\[
 J(G)\ge\int G\ge\int_{[a,x]\times[b,y]}G
 >(x-a)(y-b)\ge1/768.                                   \tag{11}
\]
Therefore a family with total cost
\(B=\sum_jw_jJ(G_j)<\infty\) has total weight of all columns that can
succeed on \(E\) at most
\[
 \sum_{j:\,E_j\cap E\ne\varnothing}w_j\le768B.         \tag{12}
\]

## 3. A fixed corner sees only a square-root slice

Fix \(\delta\in(0,\delta_0)\), let \((x,y)\in\Gamma_\delta\), and suppose
the column anchored at \((a,b)\) succeeds there.  By the preceding section
\(a,b\in(0,1)\), and therefore
\[
 2K_a(x)K_b(y)>1.                                        \tag{13}
\]
Put \(M=2u(x)u(y)\).  Since \(y\ge33/32\), direct differentiation gives
\(|u'(y)|\le4\); since \(u(x)\le3/4\),
\[
 0<M-1
 =2u(x)[u(y)-u(Y(x))]\le6\delta.                        \tag{14}
\]
Using \(K_b(y)\le u(y)\), (13) and (6) yield
\[
 u(x)-K_a(x)
 <u(x)-{1\over2u(y)}
 ={M-1\over2u(y)}\le{9\over2}\delta.                  \tag{15}
\]
Combining (3), (15), \(a>0\), and \(2x\le13/6\), gives
\[
 |a-u(x)|^2<10\delta.                                   \tag{16}
\]
Finally \(|u'(x)|\ge8/5\) on \(I\).  For fixed \(a\), (16) can hold for
an \(x\)-set of length at most
\[
 {2\sqrt{10\delta}\over8/5}
 ={5\sqrt{10}\over4}\sqrt\delta<4\sqrt\delta.         \tag{17}
\]

Let \(W_\delta(x)\) be the weighted coverage of \(\Gamma_\delta\) at its
parameter \(x\).  By (12) and (17), Tonelli gives the uniform section bound
\[
 \int_IW_\delta(x)\,dx\le3072B\sqrt\delta.             \tag{18}
\]

## 4. Consequences for exact and approximate covers

If the columns covered \(E\) almost everywhere, the unit-Jacobian
parameterization in (7) and Fubini would give
\[
 \int_IW_\delta(x)\,dx\ge|I|=1/24
\quad\hbox{for a.e. }\delta\in(0,\delta_0),             \tag{19}
\]
contradicting (18) as \(\delta\downarrow0\).  Hence an exact a.e. cover
requires infinite total column cost.

For an arbitrary finite-cost family, define its deficit
\[
 D=\int_E(1-W)_+,
\]
where \(W\) is its weighted success coverage.  From (18),
\[
 \int_I(1-W_\delta)_+\,dx
 \ge |I|-3072B\sqrt\delta.                              \tag{20}
\]
Integrating (20) over
\[
 0<\delta<h:=\min\left\{{1\over96},
                  \left({|I|\over3072B}\right)^2\right\}
\]
gives
\[
 \boxed{\displaystyle
 D\ge{|I|\over3}
 \min\left\{{1\over96},
                  \left({|I|\over3072B}\right)^2\right\}.} \tag{21}
\]
For \(B=0\), interpret the second term as \(+\infty\), in which case the
same bound is immediate.  In particular, coverage deficit cannot tend to
zero along a family with a uniform finite budget bound.

The result only rejects the strong requirement of an exact/vanishing-deficit
fixed-corner column cover of this full halo subset.  A valid strong-maximal
proof may use a different, non-columnwise budget mechanism, so (21) makes
no assertion about the optimal endpoint constant.

## 5. Compact strict-patch corollary

The obstruction is not caused by asking to cover one noncompact or boundary
receiver set.  For every \(0<\varepsilon<\delta_0/2\), let
\[
 K_\varepsilon=
 \{(x,Y(x)-\delta):x\in I,\ \varepsilon\le\delta\le\delta_0/2\}.
\tag{22}
\]
This is a compact set, and (8) makes it a compact subset of the *strict*
halo.  If a weighted fixed-corner column family covers
\(K_\varepsilon\) almost everywhere, then (18) and Fubini hold for almost
every \(\delta\in[\varepsilon,\delta_0/2]\).  Letting such depths decrease
to \(\varepsilon\) gives
\[
 B\ge {|I|\over3072\sqrt\varepsilon}
   ={1\over73728\sqrt\varepsilon}.                     \tag{23}
\]
Consequently there cannot be one budget bound independent of the compact
strict receiver set that supplies an exact fixed-corner column cover for
every such set.

This is the precise failure of the overstrong global-budget reading of the
r078 Section-6 column route: compactness of an individual receiver set does
not make an \(E\)-independent column budget available.  It does not reject
the finite fixed-cell allocation model, a receiver-local residual estimate,
or any argument not requiring a single global fixed-corner column cover of
the chosen compact set.
