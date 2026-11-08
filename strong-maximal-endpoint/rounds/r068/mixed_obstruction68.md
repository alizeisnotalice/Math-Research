# Fixed-block unbounded mixed conditional response

## Result

For every fixed \(b\ge4\), \(c>0\), and every non-node
\(s>c\), the mixed conditional multiplier derivative is unbounded:
\[
\sup_{\sum_i\xi_i=0}
\left|\partial_s\partial_c\Phi_{b,c}(s;\xi)\right|=\infty. \tag{1}
\]
Thus \(\partial_s\partial_cK_c^s\) is not a bounded \(L^2(H_b)\)
multiplier operator, even with \(b\) fixed.  This is stronger than failure
of a dimension-free \(O(b^{-1/2})\) estimate.

The proof uses only the exact normalized \(c\)-flux identity.  It does not
contradict the separate bounded first derivatives
\(\|\partial_cK_c^s\|_{2\to2}=O(1)\) and
\(\|\partial_sK_c^s\|_{2\to2}=O(b^{-1})\).

## High-frequency face

Let
\[
\xi(q)=q(b-1,-1,\ldots,-1),\qquad q>0.
\]
This is zero-sum.  Write the normalized \(c\)-flux as
\[
\partial_c\Phi=A_{b,c}(s)\{\Phi-\operatorname{Face}_c\},
\qquad
A_{b,c}(s)=bp_c\,\frac{f_{b-1,c}(s-c)}{f_{b,c}(s)}.        \tag{2}
\]
The face which removes the first, high-frequency coordinate is exactly
\[
e^{iq(b-1)c}\,
\Phi_{b-1,c}(s-c;-q,\ldots,-q)
=e^{iq(bc-s)}=:E_q(s).                                    \tag{3}
\]
There is no high-frequency decay in this term: after the first coordinate
is fixed to \(c\), all remaining frequencies are equal and their sum is
the fixed total \(s-c\).

For the full multiplier,
\[
\Phi_{b,c}(s;\xi(q))
=e^{-iqs}\,{\mathbb E}\!\left[e^{iqbZ_1}\mid S=s\right]
=O(q^{-1}).                                                \tag{4}
\]
Indeed the conditional marginal density of \(Z_1\) is
\[
z\longmapsto
\frac{q_c(z)f_{b-1,c}(s-z)}{f_{b,c}(s)}{\bf1}_{[0,s]}(z),  \tag{5}
\]
which is compactly supported and of bounded variation, away from the
finitely many knot locations.  The same bounded-variation argument gives
\(O(q^{-1})\) for every face which removes a low-frequency coordinate.
Consequently
\[
\operatorname{Face}_c={1\over b}E_q(s)+O(q^{-1}),\qquad
\partial_c\Phi=-{A_{b,c}(s)\over b}E_q(s)+O(q^{-1}).       \tag{6}
\]

## Mixed derivative

Differentiate (2) in \(s\).  The fixed-\(b\) remainder estimate needed
here is elementary and does not invoke a large-\(b\) theorem.  For any
\(n\ge3\), write the conditional one-coordinate marginal at total \(t\)
as
\[
g_{n,t}(z)=
\frac{q_c(z)f_{n-1,c}(t-z)}{f_{n,c}(t)}
{\bf1}_{[0,t]}(z).                                       \tag{7}
\]
On a compact \(t\)-neighborhood avoiding the nodes \(kc\), both
\(g_{n,t}\) and \(\partial_tg_{n,t}\) have uniformly bounded variation.
Moreover \(g_{n,t}(t)=0\), because \(f_{n-1,c}(0)=0\) for \(n\ge3\).

If a frequency vector has one coordinate \(T\) and all others \(U\), then
\[
\Phi_{n,c}(t;T,U,\ldots,U)
=e^{iUt}\int_0^t e^{i(T-U)z}g_{n,t}(z)\,dz.              \tag{8}
\]
Fourier integration by parts gives the integral \(O(|T-U|^{-1})\).
After differentiating in \(t\), the upper endpoint contributes nothing,
and the integral containing \(\partial_tg_{n,t}\) is again
\(O(|T-U|^{-1})\).  Hence
\[
\partial_t\Phi_{n,c}(t;T,U,\ldots,U)
=iU\Phi_{n,c}(t;T,U,\ldots,U)+O(|T-U|^{-1}).              \tag{9}
\]

Apply (8)--(9) first with \(n=b,T=q(b-1),U=-q\), and then to every
low-frequency face, where \(n=b-1,T=q(b-1),U=-q\) and \(t=s-c\).
This directly proves
\[
\partial_s\Phi=O(1),\qquad
\partial_s\!\left(\operatorname{Face}_c-\frac1bE_q\right)=O(1). \tag{10}
\]
The low-face frequency sum is not zero; formula (8), rather than a
sum-zero multiplier estimate, is the required argument.
The coefficient \(A_{b,c}(s)\) is positive and smooth on every interval
avoiding the finitely many nodes \(s=kc\).  Since
\[
\partial_sE_q(s)=-iqE_q(s),
\]
(2), (6), and (10) yield
\[
\partial_s\partial_c\Phi_{b,c}(s;\xi(q))
=i q\,{A_{b,c}(s)\over b}E_q(s)+O(1).                     \tag{11}
\]
The leading coefficient is strictly positive for \(s>c\).  Letting
\(q\to\infty\) proves (1). The multiplier is continuous in frequency at each finite q, so arbitrarily large values hold on sets of positive frequency measure; its essential supremum is also infinite.


## Consequence and scope

The obstruction is a moving high-frequency face phase, not a failure of
the single-\(s\) conditional \(C^1\) bounds.  It rules out the proposed
mixed-derivative route to a lossless simultaneous \((c,s)\) maximal
estimate.  It does not show that every possible mixed maximal estimate is
false, and it does not make any claim about a general strong-maximal
constant.
