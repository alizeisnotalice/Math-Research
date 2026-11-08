# r075: exact re-anchoring, and why it is not a global transport

Write `r_a=e^{-c_a}`, `r_r=e^{-c_r}`.  All rectangles below are products
coordinatewise and the receiver is in the positive cone.

## Pointwise algebra

For an actual receiver `x`, put

\[
 k={r_r-r_a\over1+r_r},\qquad
 \lambda=1-k={1+r_a\over1+r_r},\qquad a=kx,\qquad X=x-a=\lambda x. \tag{1}
\]

Then, exactly,

\[
 R_{c_a}(x)=[-r_ax,x]
 =[a-r_r(x-a),x]=a+R_{c_r}(X).                             \tag{2}
\]

There are two equivalent useful source coordinates.  If `Y=y-a`, then
`Y in R_{c_r}(X)`.  If `z=Y/lambda`, then

\[
 y=a+\lambda z,\qquad z\in R_{c_r}(x).                   \tag{3}
\]

The volume identity is also exact:

\[
 |R_{c_a}(x)|=(1+r_a)^N\prod_i x_i
 =(1+r_r)^N\prod_i X_i=|R_{c_r}(X)|.                      \tag{4}
\]

Thus, for a fixed receiver, normalized averages have no loss:

\[
 {1\over|R_{c_a}(x)|}\int_{R_{c_a}(x)}f(y)\,dy
 ={1\over|R_{c_r}(X)|}\int_{R_{c_r}(X)}f(Y+a)\,dY.       \tag{5}
\]

The same formula holds with `f` replaced by `Phi_N(f)`.  This is the exact
local compensation: the factor `lambda^N` from `dy=lambda^N dz` in (3) is
precisely the factor between the two rectangle volumes.

## Constant actual action: what is and is not transported

When `c_a` is constant, so are `k,lambda`.  On pairs set

\[
 (x,y)\longmapsto (X,Y)=(\lambda x,\ y-kx).              \tag{6}
\]

It maps the actual incidence relation `y in R_{c_a}(x)` bijectively to the
reference relation `Y in R_{c_r}(X)`.  Its Jacobian is `lambda^N`; at fixed
`x`, `dY=dy`, and (4) proves the exact kernel identity

\[
 {\mathbf1_{R_{c_a}(x)}(y)\over |R_{c_a}(x)|}\,dy
 ={\mathbf1_{R_{c_r}(X)}(Y)\over |R_{c_r}(X)|}\,dY.       \tag{7}
\]

For any nonnegative pair integrand `G`, the complete change of variables is

\[
 \int\!\int G(x,y){\mathbf1_{R_{c_a}(x)}(y)\over|R_{c_a}(x)|}
 \,dy\,dx
 =\lambda^{-N}\!\int\!\int
 G(\lambda^{-1}X,Y+(k/\lambda)X)
 {\mathbf1_{R_{c_r}(X)}(Y)\over|R_{c_r}(X)|}\,dY\,dX.    \tag{8}
\]

So the pair transformation is exact, including its receiver-Jacobian.  In
particular, it cannot be simplified by silently declaring a `lambda^N`
Jacobian in the wrong direction; for multiplicative Haar receiver measure
the displayed remaining factor disappears, whereas for Lebesgue receiver
measure it remains.

More importantly, even in this constant case (5) involves the pullback

\[
 f_X^{\sharp}(Y)=f\bigl(Y+(k/\lambda)X\bigr).            \tag{9}
\]

This is a receiver-dependent translated source, not one fixed reference
source `g(Y)`.  Consequently (7)--(8) by themselves do not transfer a
global source modular, a fixed-reference maximal theorem, or a common-anchor
statement.  They only provide an exact incidence/pair reparametrization.

## Variable action: Jacobian and pullback obstruction

For measurable `c_a(x)`, equations (1)--(5) remain valid separately at each
receiver.  There is no valid global change of variables unless the map

\[
 \Phi(x)=X=\lambda(x)x                                      \tag{10}
\]

has controlled multiplicity.  If `c_a` is `C^1`, the matrix determinant
lemma gives

\[
 D\Phi=\lambda I+x(\nabla\lambda)^T,\qquad
 \boxed{\det D\Phi=\lambda^N\bigl(1+x\cdot\nabla\log\lambda\bigr).} \tag{11}
\]

The full pair map `(x,y)->(X,Y=y-k(x)x)` is block triangular, hence has the
same determinant.  Its source pullback is now

\[
 f_x^{\sharp}(Y)=f\bigl(Y+k(x)x\bigr),                    \tag{12}
\]

which has additional receiver dependence through both the translation and
the variable dilation implicit in `X=lambda(x)x`.  For merely measurable
actions there need not even be an a.e. differentiable, one-to-one map to
which a change-of-variables theorem applies.

Critical-window amplitude control does not control the term in parentheses
in (11).  Indeed, with `r_r` fixed and `rho=sqrt(log(N)/N)`, choose on a
smooth cone-annulus

\[
 \lambda(x)=1+\varepsilon\sin(M\log|x|),\quad
 \varepsilon\asymp\rho,\quad M=\varepsilon^{-2},\qquad
 r_a(x)=(1+r_r)\lambda(x)-1.                              \tag{13}
\]

After a harmless smooth extension, `c_a(x)=-log r_a(x)` remains throughout
a critical `O(rho)` action window.  Yet

\[
 1+x\cdot\nabla\log\lambda
 =1+{\varepsilon M\cos(M\log|x|)\over1+\varepsilon\sin(M\log|x|)} \tag{14}
\]

changes sign and vanishes.  The map folds along rays.  On a logarithmic
radial interval of length comparable to `epsilon`, its output oscillation
has size comparable to `epsilon` while one period advances only
`O(epsilon^2)`; an elementary crossing count gives multiplicity
`gtrsim epsilon^{-1}` on a subinterval of output radii.  Thus a smooth,
critical-window action policy can have unbounded multiplicity
`asymp 1/rho`; no uniform global receiver-Jacobian bound follows from the
action range.

## What extra hypotheses would make this usable?

A sufficient geometric condition is radial monotonicity with quantitative
control

\[
 0<\eta\le1+x\cdot\nabla\log\lambda(x)\le M<\infty.   \tag{15}
\]

Together with `lambda` bounded above and below, this makes (10) injective
along every ray; on a suitable star-shaped cone domain it is a
ray-preserving injective `C^1` map and the area formula has the usual weight
`|det D Phi|^{-1}`.  It does **not** by itself make the map quantitatively
bi-Lipschitz on the full cone: tangential derivatives of `lambda` remain
uncontrolled.  Uniform bi-Lipschitz bounds require, in addition, a bound on
the tangential logarithmic gradient (for example a compact `C^1` bound on
`r grad_{sphere} log lambda`).  A small radial logarithmic-Lipschitz
condition, for example `|x dot grad c_a(x)|<<1`, implies (15), since

\[
 x\cdot\nabla\log\lambda
 =-{r_a\over1+r_a}\,x\cdot\nabla c_a.                    \tag{16}
\]

The required condition is positivity of the image-radius slope in (15), not
mere monotonicity of `c_a` as a radial scalar: the latter can still allow a
collapse of `r mapsto r lambda(r)`.  Even with (15), the common-source
problem (12) remains.
For a sharp global upper bound one would additionally need a source-covariant
anchor rule, or a packing/Carleson estimate that pays for the family
`f_x^sharp`.  Such hypotheses are not supplied by a measurable critical
action policy, nor by the present local-efficiency information.

## Conclusion

Re-anchoring is an exact and useful *one-receiver* identity.  For a constant
action it has an exact pair-level compensation (8), but it already produces
a receiver-indexed source pullback.  For variable actions the radial
Jacobian can fold even inside the critical window.  Therefore this algebra
does not presently yield a new sharp upper theorem; any use beyond local
identities must explicitly impose monotonicity/Jacobian control and resolve
the common-source pullback.
