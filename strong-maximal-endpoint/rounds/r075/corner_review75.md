# r075 review: tensor-Hardy charge for the fixed negative-corner policy

## Verdict

The proposed charge is correct.  For the affine policy whose actual
rectangles have the common negative corner `(-r_a,...,-r_a)`, an elementary
translation turns the whole successful halo into a one-sided tensor Hardy
halo.  Its global source-budget cost has exponential coefficient `1`,
strictly below the reference coefficient `mu>1`.  This repairs the bare
re-anchoring failure for this particular geometry; it does not supply a
general anchor-stopping theorem.

## Reduction to one-sided Hardy averages

Suppose `0<=F<=2` and the actual policy produces coordinate rectangles

\[
 [-r_a,x_i],\qquad r_a=e^{-c_a}.
\]

Put `X_i=x_i+r_a` and translate the source by

\[
 g(z)=F(z-r_a\mathbf1),\qquad z\in(0,\infty)^N.         \tag{1}
\]

Then the actual rectangle average is exactly

\[
 {1\over\prod_iX_i}\int_{[0,X]}g(z)\,dz
 =H^{\otimes N}g(X),                                      \tag{2}
\]

where `Hh(X)=X^{-1}int_0^Xh`.  Translation preserves the modular integral
`J=integral Phi_N(F)=integral Phi_N(g)`.  Thus the physical successful halo
for this policy is contained in `{H^{\otimes N}g>1}`.

This does not require an assumption that all original receiver coordinates
be positive.  The fixed-corner change of variables only requires positive
side lengths `X_i`; for an upper bound, any policy domain maps into a subset
of the full positive `X` orthant, on which (2) and the Hardy estimate apply.

## Orlicz-to-power comparison

Let

\[
 \zeta={\log(1+\log2)\over\log2},\qquad
 p=p_N=1+\zeta(N-1),\qquad p'={p\over p-1}.               \tag{3}
\]

For `0<=t<=1`, `t^p<=t=Phi_N(t)`.  For `1<=t<=2`, put
`u=log t`.  The concave function

\[
 \log(1+u)-\zeta u
\]

vanishes at `u=0,log2`, hence is nonnegative between them.  Equivalently,
`e^{zeta u}<=1+u`.  Therefore, for the usual
`Phi_N(t)=t(1+log^+t)^{N-1}`,

\[
 t^{p_N}\le\Phi_N(t)\qquad(0\le t\le2).                \tag{4}
\]

Consequently `||g||_p^p<=J`.

## Sharp enough global charge

The one-sided Hardy operator has `L^p` norm `p'`.  Tensorizing and applying
Chebyshev at level one yields

\[
 \bigl|\{H^{\otimes N}g>1\}\bigr|
 \le\|H^{\otimes N}g\|_p^p
 \le(p')^{Np}J.                                           \tag{5}
\]

The exponent has the claimed expansion:

\[
 \begin{aligned}
 Np\log p'
 &=N\bigl(1+\zeta(N-1)\bigr)
   \log\!\left(1+{1\over\zeta(N-1)}\right)\\
 &=N+{1\over2\zeta}+O(N^{-1}).                           \tag{6}
 \end{aligned}
\]

Hence this corner halo costs at most

\[
 \exp\{N+(2\zeta)^{-1}+O(N^{-1})\}\,J.                 \tag{7}
\]

Since the reference scale has logarithm `mu N-o(N)` and `mu>1`, (7) is
exponentially below that scale.  A union of `M` fixed corners costs at most
`M` times (7), so every `M=exp(o(N))` remains lower order.

## Scope

This is standard tensor-Hardy charging, not a new Hardy constant and not a
replacement for the fixed-anchor machinery.  It succeeds because all
rectangles in the policy share one physical corner.  General re-anchoring
with receiver-dependent anchors still lacks this common translation and
cannot be discharged by (2).  In particular, the result complements rather
than contradicts the affine-pushforward obstruction in `efficient75.md`.
