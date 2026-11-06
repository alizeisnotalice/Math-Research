# r055: what sharing can and cannot be forced from actual inputs

This note concerns the `sup--integral` loss in the actual strong maximal
geometry. It does not identify `C_N` with either the fixed-action model
`Q_N` or the adaptive relaxation `R_N`.

## 1. The curvature-free exact object is common-action regret

Fix a point and an outer interval `J`, and write `F_y(theta)` for the
first-coordinate average, with `theta` parametrizing an allowed interval.
Whenever the pointwise suprema admit measurable versions,
\[
 A(J)-S(J)=\inf_\theta {1\over|J|}\int_J
 \left\{\sup_\vartheta F_y(\vartheta)-F_y(\theta)\right\}\,dy.
\tag{1}
\]
This is an identity: subtract `sup_theta E F_y(theta)` from
`E sup_theta F_y(theta)`. It is the exact common-action regret and requires
neither a Hessian nor a choice of pointwise maximizing actions.

Consequently, the only curvature-free replacement for a variance hypothesis
is a direct lower bound on the right side of (1), uniformly over the outer
intervals that are near-optimal for `A`. Variance of a selected maximizer is
not invariant under ties, so it cannot by itself be a geometric theorem.

## 2. Arbitrary-dimensional flat-action counterexample

For every `n>=2`, let
\[
 f_A(x)=A\,\mathbf1_{[0,1]^n}(x),\qquad x^0=(1/2,\ldots,1/2),\qquad A>1.
\tag{2}
\]
At `x^0`, in any coordinate every interval contained in `[0,1]` and
containing `1/2` has average `A`. In the two-level notation, take an outer
interval inside the cube. Then
\[
 A(J)=S(J)=A.\tag{3}
\]
Yet choose two distinct allowable first-coordinate intervals `I_0,I_1` and
select `theta_y=I_0` on half of `J` and `theta_y=I_1` on the other half.
The selected maximizing actions have a fixed positive variance while the
regret (1) is exactly zero. Moreover the sharpness inequality
\[
 F_y(\theta_y)-F_y(\theta)\ge \kappa\|\theta-\theta_y\|^2/2
\]
fails for every positive `kappa`, since the left side vanishes for another
distinct interval inside the cube.

Thus no assertion quantified over all actual nonnegative inputs can combine
uniform fiber sharpness with a variance lower bound for "the" fiber
maximizer. Smooth plateaux give the same obstruction with arbitrarily small
curvature, so this is not an artefact of discontinuities.

More generally, for a product input `f=prod_i g_i`, actual rectangle
factorization gives exactly
\[
 M_Sf(x)=\prod_i Mg_i(x_i)=M_n\cdots M_1f(x).\tag{4}
\]
At a two-level step, `F_y(theta)=C(y)F(theta)`. All nonzero fibers have the
same set of optimizing actions; with a unique optimizer their variance is
zero, while with a tie it can be assigned arbitrarily. This is an actual
all-dimensional class on which the composition relaxation has no sharing
loss at all.

## 3. Outer near-maximizers can evade a genuine fiber mismatch

There is a second, independent obstruction to condition (3) in
`shared_action_interchange.md`: even if a larger outer interval has a strict
common-action loss, a one-sided outer interval can be equally adaptive-optimal
and have zero loss.

Here is an exact actual two-dimensional example. At receiving point
`(x_1,x_2)=(0,0)`, set
\[
 g_-=\mathbf1_{[1,2]},\qquad g_+=\mathbf1_{[2,4]},
\]
\[
 f(x,y)=g_-(x)\mathbf1_{[-1,0]}(y)+g_+(x)\mathbf1_{[0,1]}(y).
\tag{5}
\]
For an interval from `0` to `t>0`, the two first-coordinate averages are
\[
 q_-(t)=
 \begin{cases}0&t\le1,\\(t-1)/t&1<t\le2,\\1/t&t>2,\end{cases}
\quad
 q_+(t)=
 \begin{cases}0&t\le2,\\(t-2)/t&2<t\le4,\\2/t&t>4.\end{cases}
\tag{6}
\]
Both side intervals `[-1,0]` and `[0,1]` have `A=S=1/2`, using respectively
`t=2` and `t=4`. On the combined interval `[-1,1]`, adaptive choice is still
`A=1/2`, but a common interval has
\[
 S=\sup_{t>0}{q_-(t)+q_+(t)\over2}={3\over8}.\tag{7}
\]
Thus the combined outer interval has a strict mismatch, with individually
unique maximizing actions, but it cannot satisfy a uniform all-near-maximizer
loss condition: either side interval is itself exactly outer-optimal with
zero regret. Tensoring (5) with arbitrary nonnegative one-dimensional
factors gives the same physical escape mechanism in every dimension.

## 4. Endpoint-scale calibration

The flat cube and elementary product classes do not by themselves refute a
condition restricted to inputs that nearly achieve the *feedback* gap. They
are nevertheless enough to rule out an unconditional fiber theorem.

For example, the outward corner branch of the product halo of (2), restricted to every coordinate `x_i>=1`, with `L=log A`, has
volume
\[
 V_{n,+}(L)=\int_0^L {e^u u^{n-1}\over(n-1)!}\,du
 =e^L\sum_{j=0}^{n-1}{(-1)^{n-1-j}L^j\over j!}-(-1)^{n-1}.\tag{8}
\]
After division by its actual support cost `A(1+L)^(n-1)`, this branch tends
to `1/(n-1)!` as `L` tends to infinity. Even the deliberately generous
upper comparison obtained by multiplying this one branch by all `2^n`
orthants remains `exp(-n log n+O(n))`, far below the
`exp(mu n+O(sqrt(n log n)))` scale of `Q_n`. This is a calibration of the
elementary cube branch, not an upper bound for its entire strong halo or for
all product inputs.

Likewise, r044 proves that tensor products of its finite certified
two-dimensional witnesses have superexponentially decaying normalized
product-receiving score. Hence the most immediate all-dimensional
flat/product counterexamples are not known to be feedback-dangerous; they
cannot by themselves rule out a hypothesis restricted to the missing
`g log N` scale.

There is still a useful near-`Q` warning. The r028 physical fixed-ratio
construction gives actual inputs with score arbitrarily close to the fixed
action value `Q_{n,r}`. Finite disjoint translations of any compactly
truncated such witness multiply both its Orlicz cost and a disjoint union of
its valid strong receiving sets by the same factor. Thus near-`Q` score alone
cannot force a *global* action-dispersion conclusion: it is stable under
localization into arbitrarily many independent copies. This does not prove
that those witnesses have zero adaptive regret at their receiving points, so
it is not a counterexample to a condition explicitly restricted to inputs
exhibiting feedback-scale adaptive advantage.

## Conclusion

A viable rectangle-sharing theorem must be formulated in terms of the
tie-invariant regret (1), and must impose additional structural localization
excluding both flat fiber branches and side-outer-interval escapes. Neither
ordinary near-`Q` endpoint score nor a variance of an arbitrarily selected
maximizer supplies that localization. The examples above identify a real
obstruction to a universal curvature/variance route, but do not settle
whether actual feedback-dangerous near-extremizers possess enough regret to
remove the `g log N` relaxation gap.

The exploration script `shared_obstacle55.wl` was actually run in Wolfram 15.0.0. The final portable `code/shared_checks55.wl` independently checks the finite closed forms and exact slab values; its results are in `shared_checks55.json`. The exploration verified the
closed form in (8), its factorial limit, the product-average factorization,
and the exact values `1/2,1/2,3/8` in (7). Those computations corroborate
the elementary analytic derivations rather than replace them.
