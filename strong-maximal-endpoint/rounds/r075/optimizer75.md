# r075 review: compact tensor source forcing selector collapse

## Verdict

The construction is algebraically valid for all sufficiently large `N`.
It realizes the re-anchoring collapse by a bounded compact physical source
whose coordinatewise critical-window maximizing action is unique.  The
construction is deliberately not scalar-cost efficient, so it does not
contradict an efficiency-qualified packing theorem or a general sharp upper
bound.

Put `r_r=e^{-c_r}`, `ell_0=1+r_r`, `l_0=-r_r`,
`rho=sqrt(log(N)/N)`, and

\[
 H=2\rho,\quad p_r={r_r\over1+r_r},\quad
 \delta={p_r\rho\over8},\quad
 A=1+{\eta\over N},\quad b={B\over N},\quad \eta=B={1\over20}. \tag{1}
\]

Define

\[
 d={b\ell_0H\over\ell_0-2H}                              \tag{2}
\]

and let `g_N` equal

\[
 \begin{cases}
 A+b(t-l_0),&t\in[l_0-H,l_0+H],\\
 A-d,&t\in(l_0+H,1-H),\\
 A+b\ell_0+b(t-1),&t\in[1-H,1+H],\\
 0,&t\notin[l_0-H,1+H].
 \end{cases}                                               \tag{3}
\]

The endpoint jumps are harmless: this is a nonnegative bounded compactly
supported source.  Also `d=O(H/N)`, so the middle value is positive for large
`N`.

## Exact one-dimensional average

For an interval `[x-ell,x]` whose two endpoints lie in the left and right
bands of (3), respectively, write

\[
 I(x,\ell)=\int_{x-\ell}^xg_N(t)\,dt,\qquad
 a(x,\ell)={I(x,\ell)\over\ell}.                         \tag{4}
\]

At `x=1, ell=ell_0`, direct integration of the three pieces gives

\[
 \begin{aligned}
 I(1,\ell_0)
 &=AH+\frac{bH^2}{2}+(A-d)(\ell_0-2H)
   +AH+b\ell_0H-\frac{bH^2}{2}\\
 &=A\ell_0+b\ell_0H-d(\ell_0-2H)=A\ell_0.               \tag{5}
 \end{aligned}
\]

For the same endpoint regime,

\[
 \partial_x I=g_N(x)-g_N(x-\ell)=b\ell,\qquad
 \partial_\ell I=g_N(x-\ell)=A+b(x-1)-b(\ell-\ell_0).  \tag{6}
\]

Together with (5), these integrate to the exact formula

\[
 \boxed{\quad a(x,\ell)=A+b(x-1)-{b(\ell-\ell_0)^2\over2\ell}.\quad} \tag{7}
\]

The final term is strictly negative away from `ell=ell_0`; hence for each
fixed `x` the unique maximizer over every interval of admissible lengths
containing `ell_0` is `ell_0`.  Moreover
`A+b(x-1)>1` on `|x-1|<=delta`, since `delta=o(1)` and
`eta=B`.

## Endpoint containment and the maximizing actions

Let

\[
 x\in[1-\delta,1+\delta],\qquad c\in[c_r-\rho,c_r+\rho],
 \qquad \ell=x(1+e^{-c}).                                \tag{8}
\]

The right endpoint `x` lies in `[1-H,1+H]` because `delta<H`.  For the left
endpoint, write `q=e^{-c}x`.  Uniformly in (8),

\[
 q\in[r_re^{-\rho}(1-\delta),\ r_re^{\rho}(1+\delta)]. \tag{9}
\]

Since `delta/rho=p_r/8` and `r_r<1`, both deviations in (9) from `r_r` are
strictly less than `2rho=H` for all sufficiently large `N`.  Thus
`-q in[l_0-H,l_0+H]`, proving the endpoint regime used in (4)--(7) on the
entire indicated receiver/action box.

The unique action at coordinate `x` for which `ell=ell_0` is

\[
 c_*(x)=-\log\!\left({\ell_0\over x}-1\right).          \tag{10}
\]

It belongs to the action window.  Indeed `c_*(1)=c_r` and on the small
receiver interval

\[
 |c_*'(x)|={\ell_0\over x(\ell_0-x)}\le {2\over p_r}
\]

for large `N`, so `|c_*(x)-c_r|<=2delta/p_r=rho/4<rho`.

## Tensor realization and re-anchoring collapse

Set

\[
 F_N(y)=\prod_{i=1}^N g_N(y_i).                            \tag{11}
\]

It is compactly supported and

\[
 \|F_N\|_\infty
 \le\left(1+{\eta+B(\ell_0+H)\over N}\right)^N
 \le\exp\{\eta+B(\ell_0+H)\}<2                           \tag{12}
\]

for large `N` (the reference action is in the fixed positive critical
compact set, so `ell_0+H<2`).  Rectangular averaging factorizes.  Hence on
the open positive receiver cube `(1-delta,1+delta)^N`, the unique maximizing
vector action is `(c_*(x_i))_i` and the maximum average is a product of the
strict one-dimensional maxima from (7), each exceeding one.

For this actual maximizing selector,

\[
 \lambda(c_*(x_i))x_i
 ={1+e^{-c_*(x_i)}\over1+r_r}x_i
 ={\ell_0\over1+r_r}=1.                                   \tag{13}
\]

Thus its re-anchoring receiver map sends the whole open cube to the single
point `(1,...,1)`.  This is a genuine bounded compact-source realization of
the geometric collapse, not an arbitrarily prescribed policy.

## Scope

Nothing here establishes that the local `j` cost is within `J_{\rm ref}
N^{o(1)}`; the profile is intentionally a high-cost local construction.
Accordingly it does not refute an efficiency-qualified action-indexed
packing result, and it is not a counterexample to the desired global sharp
upper theorem.  It does show that uniqueness, boundedness, compact support,
and pointwise critical-window maximization alone cannot provide the missing
re-anchoring Jacobian control.

The example does not force every successful selector to collapse: other nonmaximizing successful actions may exist. Its assertion concerns the unique maximizer inside the stipulated window.
