# r079: a source-grid pricing obstruction

## Statement and scope

This is a local obstruction to pricing *all subinputs* using only anchors at
the vertices of the source grid.  It is not a counterexample to a
strong-maximal inequality, to a global \(C_N\) bound, or to a critical
efficiency statement.  In fact the original source in the example has the
displayed receiver patch in its ordinary grid-anchor halo.  The issue is
that a much cheaper successful subinput needs an interior anchor.

Let \(N\geq2\), \(1<\sigma<2\), and choose
\[
 0<d<{2-\sigma\over2\sigma},\qquad
 F=2\mathbf1_{[0,1]^N},\qquad
 E_d=(1+d,1+2d)\times(0,1)^{N-1}.                       \tag{1}
\]
The only permitted source-grid anchors are \(\{0,1\}^N\) (four of them when
\(N=2\)).  Put
\[
 u={2\sigma d\over2-\sigma}<1,\qquad
 a=(1-u,0,\ldots,0),\qquad
 G=2\mathbf1_{[1-u,1]\times[0,1]^{N-1}}.                \tag{2}
\]
Then \(0\leq G\leq F\).  For \(x\in E_d\), the positively oriented
rectangle \([a,x]\) contains all of the slab in (2), and direct cancellation
of the other coordinate lengths gives
\[
 {1\over |[a,x]|}\int_{[a,x]}G
 ={2u\over u+x_1-1}
 >{2u\over u+2d}=\sigma.                                \tag{3}
\]
The strict inequality is because the upper first-coordinate endpoint of
\(E_d\) is open.  Thus the interior anchor is genuinely successful on the
whole patch, with cost
\[
 J(G)=u\Phi_N(2)=O_{N,\sigma}(d).                        \tag{4}
\]

## Every grid-anchor column has a fixed positive charge

Let \(a'\in\{0,1\}^N\), and let a genuine source column \(0\leq G'\leq F\)
be successful at \(x\in E_d\) using the oriented rectangle \([a',x]\).
Its first anchor coordinate has to be \(a'_1=0\): if \(a'_1=1\), that
rectangle meets the support of \(F\) only on the null face \(\{y_1=1\}\),
so its \(G'\)-integral is zero.  For \(a'_1=0\), write
\[
 r_i(x)=\begin{cases}x_i,&a'_i=0,\\1-x_i,&a'_i=1,
          \end{cases}\qquad (i\geq2).
\]
Success and \(\Phi_N(t)\geq t\) give
\[
 J(G')\geq\int G'\geq\int_{[a',x]}G'
 \geq \sigma x_1\prod_{i=2}^N r_i(x)
 \geq \sigma x_1\prod_{i=2}^N\min\{x_i,1-x_i\}.          \tag{5}
\]

Consequently, let \(\{(w_j,a_j,G_j)\}\) be any finite or countable legal
weighted cover of \(E_d\) by such genuine grid-anchor columns: the sum of
the weights of columns successful at each \(x\) is at least one.  At any
such \(x\), (5) implies
\[
 \sum_jw_jJ(G_j)
 \geq \sigma x_1\prod_{i=2}^N\min\{x_i,1-x_i\}.
\]
Taking \(x_1\uparrow1+2d\) and \(x_i\to1/2\) yields
\[
 \boxed{\displaystyle
 \sum_jw_jJ(G_j)\geq\sigma(1+2d)2^{-(N-1)}.}            \tag{6}
\]
If coverage is required only almost everywhere, apply (5) first on a
positive-measure neighborhood of these limiting points, choose a covered
point there, and then shrink the neighborhood.  Thus (6) is unchanged.
The lower bound is independent of \(d\), whereas (4) tends to zero.

For contrast, the fixed full source \(F\) itself is grid-anchor successful
on \(E_d\): with first anchor zero (and either choice in the other
coordinates), its average is \(2/x_1>\sigma\).  Thus this is not a claim
that the fixed-\(F\) halo needs an interior anchor.

## Exact whole-patch grid price

There is an even sharper one-column formulation.  Fix any
\(a'\in\{0,1\}^N\) with \(a'_1=0\).  If \(G'\leq F\) succeeds on all of
\(E_d\), send \(x_1\uparrow1+2d\) and, for each \(i\geq2\), send \(x_i\)
to the endpoint opposite \(a'_i\).  The limiting rectangle contains the
whole unit cube and has volume \(1+2d\).  By continuity of the rectangle
integral,
\[
 \int_{[0,1]^N}G'\geq\sigma(1+2d).                       \tag{7}
\]
Jensen on the unit cube gives
\[
 J(G')\geq\Phi_N\!\bigl(\sigma(1+2d)\bigr).              \tag{8}
\]
The bound is attained by the uniform real column
\[
 G'_0=\sigma(1+2d)\mathbf1_{[0,1]^N},
\]
which lies below \(F\) by (1) and has average
\(\sigma(1+2d)/x_1>\sigma\) on every rectangle in the patch.  Hence
\[
 c_{a'}(E_d)=\Phi_N\!\bigl(\sigma(1+2d)\bigr)
 \quad(a'_1=0),
 \qquad c_{a'}(E_d)=+\infty\quad(a'_1=1).                \tag{9}
\]
In particular, the ratio of the best whole-patch source-grid price to the
interior-anchor price (4) is at least
\[
 {\Phi_N(\sigma(1+2d))\over u\Phi_N(2)}\longrightarrow\infty
 \qquad(d\downarrow0).                                  \tag{10}
\]
This proves that source-grid vertex pricing cannot approximate every
successful subinput with a uniform multiplicative cost factor, even though
the original fixed source has an exact grid-anchor halo representation on
this patch.

## Exact continuous-anchor price for one real column

The preceding gap is not an artifact of using the particular slab in (2).
For this paragraph, \(H_aG(x)\) denotes the average over the oriented
rectangle with anchor \(a\) and receiver \(x\), whenever that rectangle is
nondegenerate.  Let \(C_{\rm one}(N,\sigma,d)\) denote the infimum of
\(J(G)\) over one continuous anchor and one \(0\leq G\leq F\) such that
\(H_aG(x)\geq\sigma\) for almost every nondegenerate \(x\in E_d\).
Equivalently here one may require \(>\sigma\): the optimizing slab below is
strict on the open patch.  If a transverse anchor coordinate lies in
\((0,1)\), its receiver hyperplane \(x_i=a_i\) is simply omitted as a
Lebesgue-null degeneracy set; the realizing optimizer can instead take all
transverse anchor coordinates equal to zero, so it has no such degeneracy.
For this one-column problem all continuous anchors and all subinputs can be
solved exactly.  First project every transverse anchor coordinate into
\([0,1]\): if it lies below zero or above one, replacing it by the nearer
cube face leaves its intersection with \([0,1]^N\) unchanged and only
shortens every receiver rectangle.  The first coordinate must be below one;
if it is negative, replace it by zero by the same argument.  Thus write it
as \(a_1=1-u\), \(0<u\leq1\).  Anchors with \(a_1\geq1\) have zero source
intersection on this receiver patch.

Fix \(u\) and transverse coordinates \(a_i\in[0,1]\).  Splitting the cube
at those transverse coordinates produces \(2^{N-1}\) orthants \(Q_\eta\)
(zero-volume pieces are discarded).  For each \(Q_\eta\), let the
corresponding receiver coordinates approach the endpoint opposite \(a_i\);
also let \(x_1\uparrow1+2d\).  The limiting rectangle intersects the source
cube in exactly
\[
 S_\eta=[1-u,1]\times Q_\eta
\]
and strict success on \(E_d\), followed by continuity, forces
\[
 \int_{S_\eta}G\geq \sigma(u+2d)|Q_\eta|.               \tag{11}
\]
Jensen on each positive-volume \(S_\eta\), followed by summation, gives
\[
 J(G)\geq
 u\sum_\eta |Q_\eta|
 \Phi_N\!\left(\sigma\left(1+{2d\over u}\right)\right)
 =u\Phi_N\!\left(\sigma\left(1+{2d\over u}\right)\right). \tag{12}
\]
The argument also handles boundary anchors \(a_i=0,1\): their zero-volume
orthants are omitted and the remaining orthants still have total volume one.

Conversely, for any \(u\) with
\[
 g:=\sigma\left(1+{2d\over u}\right)\leq2,               \tag{13}
\]
the column
\[
 G_u=g\mathbf1_{[1-u,1]\times[0,1]^{N-1}}
\]
lies below \(F\) and has rectangle average
\[
 {gu\over u+x_1-1}>{gu\over u+2d}=\sigma
\]
at every \(x\in E_d\), independently of the transverse anchor coordinates.
Thus (12) is exact.  Eliminating \(u=2\sigma d/(g-\sigma)\) proves the
continuous one-column formula
\[
 \boxed{\displaystyle
 C_{\rm one}(N,\sigma,d)
 =\min_{\sigma(1+2d)\leq g\leq2}
 {2\sigma d\,\Phi_N(g)\over g-\sigma}.}                 \tag{14}
\]
Here the lower endpoint is the condition \(u\leq1\), and the upper endpoint
is precisely the cap \(G_u\leq2\).

For completeness, on \(g>1\) the derivative of the objective in (14) has
the sign of
\[
 D(g)=(N-1)g-\sigma(N+\log g),\qquad
 D'(g)=N-1-{\sigma\over g}>0\quad(g>\sigma).             \tag{15}
\]
Hence the minimizer is the unique zero of \(D\), clipped to the interval in
(14).  When \(N=2\), \(D(2)=2-\sigma(2+\log2)<0\), so the cap endpoint
\(g=2\) is the minimizer.

One useful asymptotic makes the obstruction quantitative.  Let
\(\sigma_N=1+\kappa/N\), with fixed \(\kappa>0\), and
\(d=N^{-\beta}\), \(\beta>1\).  The interior zero of (15) has
\[
 g_*=1+{\kappa+1\over N}+o(N^{-1}),\qquad
 u_*={2\sigma_Nd\over g_*-\sigma_N}\sim2Nd,
\]
so
\[
 C_{\rm one}(N,\sigma_N,d)\sim2e^{\kappa+1}Nd.
\]
By (9), the whole-patch source-grid price is
\[
 \Phi_N\!\left(\sigma_N(1+2d)\right)\sim e^\kappa.
\]
Their ratio and logarithmic loss are therefore
\[
 {\,\Phi_N(\sigma_N(1+2d))\over C_{\rm one}(N,\sigma_N,d)}
 \sim {N^{\beta-1}\over2e},\qquad
 \log{\Phi_N(\sigma_N(1+2d))\over C_{\rm one}(N,\sigma_N,d)}
 =(\beta-1)\log N-\log(2e)+o(1).                        \tag{16}
\]
This is an exact continuous-anchor **one-column** calculation only.  It
does not solve a full fractional-capacity allocation problem and gives no
general \(C_N\) conclusion.
