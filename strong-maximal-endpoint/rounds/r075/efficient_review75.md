# Audit of the qualified nonsingular reanchoring example

## Verdict

**PASS as a local, qualified geometric obstruction, with one size
qualification.**  Starting from the compact capped r074 source at
\(x_0=(1,\ldots,1)\), the proposed smooth injective vector policy preserves
strict success and the local cost on some nonempty open receiver set.  Its
canonical reanchor has no folding but its weighted pushforward density is
\(\mu^{N+o(N)}\), so a subpolynomial vector-reanchor distortion cannot be
deduced from smoothness, injectivity, actual success, and the local-cost
screen alone.

This is a statement about the proposed weighted geometric reanchoring
interface.  It does not prove a norm inequality for \(j_W^{c(x)}\), does
not assert that the prescribed policy maximizes the actual action, and does
not make an aggregate halo claim.

## 1. Local persistence of the qualified source

Let \(F\) be the r074 compact physical source, with \(0\le F\le2\),
support away from coordinate hyperplanes, strict success at \(x_0\) under
the scalar action \(c_a\), and
\[
 j_W^{c_a}\Phi_N(F)\le J_{\rm ref}N^{o(1)}.              \tag{1}
\]
For a receiver \(x\) near \(x_0\), set coordinatewise
\[
 c_i(x)=c_a+\log x_i,\qquad r_a=e^{-c_a}.
\tag{2}
\]
Then the actual interval in coordinate \(i\) is exactly
\[
 [-e^{-c_i(x)}x_i,x_i]=[-r_a,x_i].                        \tag{3}
\]
For \(|x_i-1|\le N^{-2}\), \(|\log x_i|\le2N^{-2}\), which is
\(o(\rho)\).  Since \(c_a-c_r=O(\rho L^{-1/4})\), all these actions remain
in the critical window for large \(N\).

The received average is continuous in \(x\): write it as the normalized
integral of the bounded compact function \(F\) over the product in (3),
then use dominated convergence.  The same is true for the local cost.  One
explicit physical expression for the latter has the form
\[
 j(x)={1\over\prod_i x_i}\int_{\prod_i[-r_a,x_i]}
 {\Phi_N(F(y))\over
 H_{\mathbf c(x)}\!\left(\sum_i\log{x_i\over|y_i|}\right)}\,dy.
\tag{4}
\]
On the compact source support this integrand has a common local dominator.
Indeed the sum in (4) stays in a compact subset of \((0,\infty)\), and
the zero-jump term gives
\(H_{\mathbf c}(s)\ge s^{N-1}/\Gamma(N)>0\).  Thus DCT applies to (4) as
well.

It follows that, for every fixed \(\varepsilon>0\), there is an open
nonempty \(E_N\) containing \(x_0\), which may be taken inside
\([1-N^{-2},1+N^{-2}]^N\), such that
\[
 \text{received average}>1,\qquad
 j(x)\le J_{\rm ref}N^\varepsilon\quad(x\in E_N).
\tag{5}
\]
The containment rather than equality with the entire displayed cube is the
strict conclusion from DCT alone.  To claim the whole \(N^{-2}\)-cube one
would need a quantitative modulus of continuity compared with the chosen
strict-success margin.  This distinction does not affect the geometry
below, which holds on every such open set.

## 2. Smooth injective reanchor and its weighted density

Let \(r_r=e^{-c_r}\), and use the coordinatewise reference reanchor
\[
 \lambda_i(x)={1+r_a/x_i\over1+r_r},
 \qquad T_i(x)=\lambda_i(x)x_i={x_i+r_a\over1+r_r}.
\tag{6}
\]
This map is affine and injective, with
\[
 \det DT=(1+r_r)^{-N}.                                   \tag{7}
\]
Thus the previous rank-collapse mechanism is absent.  If the intended
source-pair/Fubini cancellation carries the natural rectangle-scale weight
\(w(x)=\prod_i\lambda_i(x)\), its density after changing variables is
exactly
\[
 {w(x)\over|\det DT(x)|}
 =\prod_{i=1}^N\left(1+{r_a\over x_i}\right).            \tag{8}
\]
On \(E_N\), \(x_i=1+O(N^{-2})\), while
\(r_a=e^{-\mu+O(\rho)}\) and \(1+e^{-\mu}=\mu\).  Therefore
\[
 \log {w\over|\det DT|}
 =N\log\mu+O(N\rho)+O(N\cdot N^{-2})
 =N\log\mu+o(N),
\tag{9}
\]
or equivalently
\[
 {w\over|\det DT|}=\mu^{N+o(N)}.                         \tag{10}
\]
So even a smooth affine one-to-one reanchor has exponentially large
weighted distortion on a set whose receivers are actually strict and pass
the fixed-\(\varepsilon\) local-cost screen.

## Scope of the weight calculation

Equations (8)--(10) are conditional on the reanchoring scheme genuinely
requiring the weight \(w=\prod_i\lambda_i\); they are a bare geometric
pushforward calculation.  They do not establish that an arbitrary physical
source cost, or \(j_W^{c(x)}G\), has this norm or that no different
action-indexed packing can succeed.  The set \(E_N\) may be extremely
small, so it creates no aggregate contradiction.  These qualifications are
essential: the example identifies a strong obstruction to a particular
vector reanchor/Fubini route, not to a general endpoint upper bound.
