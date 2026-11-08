# Reanchoring variable actions: an exact Jacobian obstruction

## Conclusion

A shrinking critical action window by itself gives **no** uniform lower
Jacobian, finite-multiplicity, or unweighted-pushforward transport estimate
for the
reanchoring map
\[
 T(x)=\lambda(c(x))x,\qquad
 \lambda(c)={1+e^{-c}\over1+e^{-c_r}}.
\tag{1}
\]
This remains false for a smooth, radially monotone policy.  There is an
explicit policy in the critical window which sends a positive-measure
annulus exactly onto a sphere.  Thus a proof which first reanchors every
variable-action receiver and then applies a fixed-reference Fubini bound
needs a genuinely quantitative distortion/packing hypothesis.

This is a geometric obstruction to the bare-Jacobian/Fubini reanchoring
interface, not a counterexample to an aggregate strong-maximal upper bound
or to a finite action-envelope estimate.

## 1. Exact smooth collapse inside the critical window

Put \(\rho=\sqrt{\log N/N}\), fix \(A>0\), and restrict actions to
\[
 I_N=[c_r-A\rho,c_r+A\rho].
\tag{2}
\]
The function \(\lambda\) is smooth and strictly decreasing.  Write
\(\lambda_-=\lambda(c_r-A\rho)>1\) and
\(\lambda_+=\lambda(c_r+A\rho)<1\).  For any \(r_0>0\), let
\[
 \mathcal A_N=\{x\in(0,\infty)^N:\ r_0/\lambda_-<|x|<r_0/\lambda_+\}.
\tag{3}
\]
It has positive \(N\)-dimensional measure.  On this annulus define the
smooth radial policy
\[
 c(x)=c(|x|),\qquad
 c(r)=-\log\!\left((1+e^{-c_r}){r_0\over r}-1\right).
\tag{4}
\]
The expression in parentheses is positive on (3), and direct inversion
shows
\[
 \lambda(c(r))={r_0\over r},\qquad
 c(r_0/\lambda_-)=c_r-A\rho,\qquad
 c(r_0/\lambda_+)=c_r+A\rho.
\tag{5}
\]
Hence \(c(\cdot)\) takes values exactly in the prescribed shrinking
window.  It is also increasing:
\[
 c'(r)={(1+e^{-c_r})r_0/r^2
 \over (1+e^{-c_r})r_0/r-1}>0.
\tag{6}
\]
Nevertheless,
\[
 T(x)={r_0\over|x|}x\in r_0\mathbb S^{N-1}\cap(0,\infty)^N,
 \qquad x\in\mathcal A_N.
\tag{7}
\]
Thus a smooth monotone policy collapses the whole annulus to a
codimension-one set.  The window is small but nontrivial: its relative
radial width is \(\lambda_+^{-1}-\lambda_-^{-1}=\Theta(\rho)\).

## 2. The correct Jacobian and area formula

For an arbitrary differentiable policy,
\[
 DT=\lambda I+x\otimes\nabla\lambda,
 \qquad
 \det DT=\lambda^N\bigl(1+x\cdot\nabla\log\lambda\bigr).
\tag{8}
\]
For a radial policy this becomes, with \(t(r)=r\lambda(c(r))\),
\[
 \det DT=\lambda^{N-1}t'(r)
 =\lambda^{N-1}\{\lambda+r\lambda'(c)c'(r)\}.
\tag{9}
\]
The policy (4) has \(t(r)=r_0\), so (9) vanishes identically.  This is the
precise failure; smoothness and monotonicity of \(c\) do not address it.

For a Lipschitz map, the valid area formula is
\[
 \int_E g(Tx)|\det DT(x)|\,dx
 =\int_{\mathbb R^N}g(a)N(T,E,a)\,da,
\tag{10}
\]
where \(N(T,E,a)\) is the multiplicity.  If \(T\) is nonsingular a.e.
and has controlled multiplicity, the unweighted pushforward instead is
\[
 \int_Eg(Tx)\,dx
 =\int_{\mathbb R^N}g(a)
   \sum_{x\in T^{-1}(a)\cap E}{1\over|\det DT(x)|}\,da.
\tag{11}
\]
There is no version of (11) with a bounded density for (7): the pushforward
of Lebesgue measure on \(\mathcal A_N\) is singular on the sphere.  For
example, if \(g_\epsilon\) is a continuous cutoff equal to one on the
\(\epsilon/2\)-tube of \(r_0\mathbb S^{N-1}\) and supported in its
\(\epsilon\)-tube, then
\[
 \int_{\mathcal A_N}g_\epsilon(Tx)\,dx=|\mathcal A_N|,
 \qquad \int g_\epsilon(a)\,da=O_N(\epsilon).
\tag{12}
\]
No uniform \(L^1\) reanchoring estimate follows from the action-window
amplitude alone.

## 3. What a sufficient monotone policy condition actually is

Mere radial monotonicity is insufficient by (6)--(7).  A usable radial
condition is, for some fixed \(\delta>0\),
\[
 t'(r)=\lambda+r\lambda'(c)c'(r)\ge\delta\lambda
\tag{13}
\]
on the receiver interval, together with injectivity of \(t\).  Then
\[
 |\det DT|\ge\delta\lambda^N,
\qquad
 \int_Eg(Tx)dx
 \le\delta^{-1}\sup_E\lambda^{-N}\int_{T(E)}g(a)da.
\tag{14}
\]
For nonradial policies, (13) is replaced by
\[
 1+x\cdot\nabla\log\lambda(c(x))\ge\delta
\tag{15}
\]
plus a finite multiplicity bound.  These are genuine additional *local*
policy regularity/packing assumptions, stated on a compact receiver
annulus rather than as a global finite-distortion theorem.

Even under (13), the critical window only gives
\(\lambda^{\pm N}=\exp(O(N\rho))=\exp(O(\sqrt{N\log N}))\).  Therefore
an *unweighted* receiver-volume change is far too expensive at the desired
logarithmic precision.  One needs either an exact compensating weight or an
action window of order \(O((\log N)/N)\) to reduce that factor to a power of
\(N\).

## 4. Why rectangle-volume cancellation is not enough

The special choice (1) has the exact one-rectangle volume identity
\[
 (1+e^{-c})^N\prod_i x_i
 =(1+e^{-c_r})^N\prod_i[\lambda(c)x_i].
\tag{16}
\]
This can cancel \(\lambda^N\) in a Fubini calculation that carries the
full unnormalised rectangle-volume weight.  It does not cancel the
derivative factor in (8).  Indeed, after the formal cancellation, changing
variables still leaves
\[
 {1\over1+x\cdot\nabla\log\lambda(c(x))},
\tag{17}
\]
as well as multiplicity.  In the exact collapse (4), this denominator is
zero.  Hence volume matching alone cannot make a variable reanchoring into
a measure-preserving source pairing.

The action-indexed local cost \(j_W^{c(x)}G\) is normalized conditional
source cost.  A fixed-action Fubini/Markov estimate may use its stochastic
kernel, but inserting \(c(x)\) first creates exactly the pushforward issue
above if it is reduced to a fixed \(c_r\) source coordinate.  The known
positive action envelope still gives a finite, \(N\)-dependent bound; this
note does not challenge it.  It only says that the formal \(\lambda^N\)
cancellation cannot improve that argument by a bare unweighted change of
variables, because (12) supplies no such pushforward bound.

Nor has the tubular test \(g_\epsilon\) in (12) been realized here as a
legal physical source cost \(G\) satisfying the efficiency qualification.
It therefore does not refute an actual-policy global cost estimate.  A
valid improved global estimate may retain an action-indexed source
allocation/stopping mechanism, or prove a distortion condition from
additional structure; it simply cannot obtain the required reanchoring
step from \(|c-c_r|=O(\rho)\) alone.

## 5. Conditional additive loss ledger (not an actual-policy theorem)

This is only algebra under an additional, unproved common-source pullback
estimate.  Suppose that estimate costs
\(D_N=M_N/\delta_N=N^{\beta+o(1)}\), multiplying both the good captured
band transport and the bad wide-window screen.  With screen threshold
\(J_{\rm ref}N^\varepsilon\), their relative powers are
\[
 \ell_{\rm good}=K-\kappa+\beta,
 \qquad \ell_{\rm bad}=\tfrac12+\beta-\varepsilon,
 \qquad \varepsilon<2K,\quad K>\kappa.
\tag{18}
\]
The two contributions are added.  Take \(K\downarrow\kappa\) and choose
\(\varepsilon\in(1/2,2K)\) sufficiently close to \(2K\).  Then
\[
 \ell_{\rm good}=\beta+o(1),\qquad
 \ell_{\rm bad}=\beta+\tfrac12-2\kappa+o(1)<\beta+o(1),
\tag{19}
\]
because the shifted-reference \(\kappa=\alpha=(9-\mu)/24>1/4\).
Thus this additive accounting has total relative loss \(N^{\beta+o(1)}\).
It recovers coefficient \(\kappa\) only for subpolynomial \(D_N\).
Neither this ledger nor the geometric slope conditions prove a general upper
bound: the common-source pullback estimate remains the unresolved
interface.
