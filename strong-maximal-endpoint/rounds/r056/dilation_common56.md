# r056 — common-parameter dilation identity

## Scope

This is the full-positive-quadrant, one-receiver, two-equal-block, regular
KKT calculation with cap two and target one.  It is a common path
`c_1=c_2=c`, not the contrast path.  It makes no statement about multiple
receivers, finite moving boxes, or `C_N`.

## Fixed-shape pullback

Put `N=2b`, `t_i=cx_i`, and `S=x_1+x_2`.  The block likelihood has the form
\[
R_{b,c}(cx)={c^{b-1}\over\Gamma(b)Z(c)^b}S_b(x),
\]
so the `x`-shape is fixed and
\[
\rho:=(\log R)'={Nm(c)-2\over c}.
\]
It is generally nonzero.  The pulled-back base and likelihood measures have
scores
\[
L_\xi={2\over c}-S,\qquad L_\nu={Nm(c)\over c}-S=L_\xi+\rho,
\]
and logarithmic second derivatives
\[
L_\xi'=-{2\over c^2},\qquad L_\nu'=-{Nv(c)\over c^2}.
\]
All four algebraic formulas were checked by the local Wolfram run in
`dilation_common56.wl`.

Let `a=lambda R`, `H=h(a)`, and define
\[
A=E_\nu[a h'(a)]>0,\qquad B=E_\nu[H L_\nu].
\]
The strict inequality is the full regular-response qualification: it includes
the soft derivative and the lower-coarea contribution at the KKT jump.  If
\(\eta=(\log a)'\), differentiation of `E_nu H=1` gives
\[
\eta=-{B\over A}.
\]
Thus a common change cannot be removed merely by observing that `R` has a
fixed shape: the multiplier changes through this scalar response.

## Exact first and second derivatives

Write \(q(a)=a-\Phi^*_{[0,2]}(a)\).  Fenchel equality and target one give
`P=int q dxi` and `q'(a)=1-H`.  Since
\(\int a(1-H)d\xi=\lambda E_\nu(1-H)=0\),
\[
P'=\int qL_\xi\,d\xi.
\]
For the second derivative, use
\[
\int a(1-H)L_\xi d\xi
=\lambda E_\nu[(1-H)(L_\nu-\rho)]= -\lambda B.
\]
This gives the exact common formula
\[
\boxed{\quad
P''=\int q\left(L_\xi^2-{2\over c^2}\right)d\xi
       +{\lambda B^2\over A}.
\quad} \tag{1}
\]
The last term is a nonnegative Schur response.  In `t` coordinates, set
\(T=t_1+t_2=cS\) and \(C_T=\operatorname{Cov}_\nu(H,T)=-cB\).  Then
\[
\boxed{\quad
c^2(\log P)''=
 {\int(T-2)^2q\,d\xi\over P}-2
-\left({\int(T-2)q\,d\xi\over P}\right)^2
+{\lambda\over P}{C_T^2\over A}.
\quad} \tag{2}
\]
Equivalently, the first three terms are the signed-`q/P` variance of `T`
minus two.  The formula remains valid away from a common stationary point;
the displayed square is the usual `-(P'/P)^2` term.

Because `nu` is a product measure, `H` is coordinatewise nondecreasing and
`T` is increasing, FKG gives `C_T>=0`.  This confirms the response sign but
does not control its size.

## What this does not yield

The r054 contrast theorem alone controls a signed **difference** moment of
`t_1-t_2`; formula (2) alone does not determine the radial signed variance
of `T` or the Schur quotient `lambda C_T^2/(PA)`.  Thus away from a common
stationary point the two `M`-dependent terms in (2) remain uncontrolled by
this reduction.

At a common stationary point, however, the companion common-path analysis
uses the original-parameter Hessian together with the r054 `Z`-plus closure,
the r053 response estimate, and Cauchy--Schur control to establish the
leading `Np_mu` curvature.  The present note neither reproves nor weakens
that stationary theorem: it records the full nonstationary moment formula
that explains why its additional stationary input is needed.  There is still
no conclusion here for multiple receivers or `C_N`.

## Local computation record

Wolfram 15.0.0 for Mac OS X ARM returned `True` for the derivative of the
common likelihood scalar, the score difference, the `nu` second derivative,
and the formal Schur cancellation.  The Wolfram computation corroborates
algebra only; (1)--(2) follow by the displayed differentiation.
