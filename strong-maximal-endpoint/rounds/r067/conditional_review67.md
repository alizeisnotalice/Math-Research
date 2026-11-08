# Independent audit of conditional_flux67

## Verdict: PASS, with two presentation clarifications

The conditional labelled first-response theorem is mathematically
well-founded in its stated near-\(\mu\), critical-sum scope.  As an
\(L^2\) Fourier multiplier statement it needs no height cap; a fixed cap
enters only when it is paired with the separate Bregman-budget argument.
The random-interval lemma supplies the missing \(b^{-1/2}\) numerator
bound, the fixed mean-matched tilt is differentiated correctly, and the
joint sign factors do not lose the bound.  The theorem remains a
conditional angular \(L^2\) statement; it does not solve the global
source-budget or receiver-cone packing problem.

Two details should be stated explicitly if this is imported elsewhere:

1. all differentiation is at a fixed chosen tilt \(\theta\), not along
   the mean-matching curve \(c\mapsto\theta(c)\); and
2. the finite-\(V\) application must use zero extension of \(F-B\), not
   zero extension of the constant \(B\).

Neither is a gap in the current proof.

## 1. Normalized flux and conditioning

For
\[
 q_c(z)={e^{-z}(1+1_{\{z>c\}})\over1+e^{-c}}1_{\{z\ge0\}},
\]
distributional differentiation gives
\[
 \partial_cq_c=p_cq_c-p_c\delta_c.
\]
Differentiating the raw \(b\)-fold convolution and its density, then
using the quotient rule, gives exactly equation (1) of the theorem.  The
factor \(b\) in
\[
 A_{b,c}=bp_c\,{f_{b-1,c}(s-c)\over f_{b,c}(s)}
\]
comes from the density derivative and is correct.  The face-recentering
identity (2) also checks directly from
\(\sum_j\xi_j=0\).

At the selected point, \(m_{c,\theta}=s/b\), so
\[
 \prod_j e^{-im(\xi_j+\tau)}
 =e^{-is\tau}
\]
on the sum-zero frequency hyperplane.  Thus \(I_\xi/I_0\) is the
conditional characteristic function with the correct Fourier-inversion
normalization, and
\[
 I_0=2\pi f^{(\theta)}_{b,c}(s)\asymp b^{-1/2}.
\]
The density lower bound follows uniformly from the tilted exponential
factor and the standard local CLT at the tilted mean.

The apparent centering issue is handled correctly: the \(c\)-derivative
in the proof is the derivative of the raw fixed-\(\theta\) numerator and
denominator.  The factor \(e^{-imt}\) is used only to represent their
values at the evaluation point.  Consequently the raw score contains
\(i(\partial_cm)t\), which is precisely the term removed by the
\(\tau\)-integration by parts.  No derivative \(\theta'(c)\) is allowed
or needed, since conditioning on the total is invariant under a fixed
exponential tilt.

## 2. Random-interval estimate and integration by parts

The shifted-product lemma (3) is valid uniformly in all shifts.
The layer-cake identity (5) is exact: after deleting interval \(I_i\),
the integrand is the expected increment of intersection length.
At most two intervals can be extremal (unique maximal left endpoint or
unique minimal right endpoint), and every relevant leave-one-out
intersection has length at most \(2R_{(2)}\).  The survival formula (6)
is the correct event that at most one radius lies below \(r\).  Its
small-\(r\) and large-\(r\) estimates yield \(E R_{(2)}=O(b^{-1/2})\).
Together with \(h\,w\lesssim1-w\), this proves (3).

For \(b\ge3\), the product has enough decay at infinity to justify both
the differentiated Fourier integral and
\(\int P_\xi'(\tau)d\tau=0\).  Equations (7)--(8) have the right
residual sizes: Taylor expansion gives \(O(t^2)\) at zero, while the
global logarithmic derivatives grow by at most \(O(|t|)\).  Applying
(3) to (9) therefore gives
\[
 |\dot I_\xi|+|\dot I_0|=O(b^{-1/2}),
\]
and division by \(I_0\) proves the dimension-free bound (T).

## 3. Joint labels

For a sign character, the kernel
\[
 q_{c,\varepsilon}(z)=
 {e^{-z}(1+\varepsilon1_{\{z>c\}})\over1+e^{-c}}1_{\{z\ge0\}}
\]
is correctly understood as a signed all-plus expectation, not a
probability density.  Its normalization denominator must remain
all-plus, as in equation (12).  The quotient is nevertheless bounded by
one, because it is the all-plus conditional expectation of a Fourier
character times the product of the \(\pm1\) jump labels.

For minus factors, the small-frequency contraction (14) is uniform:
the action window stays in a compact subset of \((0,\infty)\), so
\[
 {|1-re^{ict}|\over1+r}\le\kappa<1
\]
when \(|t|\le\pi/(2c_{\max})\).  Stratifying by the number of
marked-small minus factors gives \(m\kappa^m\), which is uniformly
bounded; Hölder gives
\[
 \int\prod_jw_\gamma(\tau+x_j)\,d\tau=O(b^{-1/2}).
\]
For the remaining minus factors, (16) reduces their
\(O(1+|t|)\) residual to the same random-interval estimate.  This
validly gives the labelled analogues of (9)--(10) and hence (17).

## 4. Exact finite-depth use

The conditional theorem applies without a tangent-limit approximation.
For one block set
\[
 X'=V\xi-V{\bf1}/b,\qquad
 Y'=(V+s)\eta-(V+s){\bf1}/b.
\]
Since \(\eta=(V\xi+Z)/(V+s)\) conditional on \(\sum Z_i=s\),
\[
 Y'=X'+\bigl(Z-s{\bf1}/b\bigr)=X'+PZ.                    \tag{A}
\]
For a radial scalar \(B(w)\), the source constant is mapped to the same
constant by every action.  Hence the angular comparison is applied to
the zero extension of
\[
 G(\eta)=F(w,\eta)-B(w)
\]
from the source simplex to the sum-zero hyperplane; it is not applied to
a zero extension of \(B\) itself.  The output may then be restricted back
to the receiver simplex.

The uniform Dirichlet push gives the exact norm comparison
\[
 \|G\circ T_{V,Z}\|_{L^2(d\sigma_\xi)}^2
\le\left({V+s\over V}\right)^{b-1}
\|G\|_{L^2(d\sigma_\eta)}^2.                              \tag{B}
\]
Multiplication by the receiver radial Jacobian \(V^{b-1}\) changes the
right side to \((V+s)^{b-1}\), the source radial Jacobian.  Thus there is
no unrecorded finite-\(V\) Jacobian loss in this conditional \(L^2\)
comparison.  Source and receiver radial domains are still different;
this cancellation is local in \((V,s)\), not a global cone allocation.

## 5. Conditional radial sensitivity: PASS

The labelled radial distributional identity is correct:
\[
 \partial_zq_{c,\varepsilon}
 =-q_{c,\varepsilon}+a_0\delta_0+\varepsilon p_c\delta_c.
\]
For an \(s\)-derivative of a convolution one differentiates one selected
factor, or equivalently averages that identity over the \(b\) possible
factors.  This explains why the normalized face coefficients in RF2 are
\(\alpha_0,\alpha_c\), rather than \(b\alpha_0,b\alpha_c\).  The
all-plus denominator obeys the same formula, so quotienting yields RF2.
RF3 is the earlier action flux, and RF4 follows by differentiating the
Gamma factor conditional on the posterior jump count:
\[
 \partial_s\log f_{b,c}(s)
 =-1+(b-1){\mathbb E}[ (s-cK)^{-1}\mid S=s].
\]
There is therefore no omitted normalizer derivative or sign coefficient.

For \(\partial_s\), the raw fixed-tilt Fourier numerator has derivative
\[
 -i\int\tau P_\xi(\tau)\,d\tau.
\]
Equation (9) and the shifted-product estimate give
\[
 \left|\int\tau P_\xi\right|=O(b^{-3/2});
\]
division by \(I_0\asymp b^{-1/2}\), including the differentiated
all-plus denominator, yields
\[
 \|\partial_s\mathcal K_c^s\|_{2\to2}=O(b^{-1}).
\]
The marked-small/marked-large sign proof supplies the same bound for
every joint label character.  The mean-matched tilt remains only a
pointwise representation of the raw integrand.  No \(\theta'(s)\) term is
present or permissible.

The path estimate in equation (25) is a Hilbert-space statement for a
fixed function.  Equation (26) correctly applies both kernels to the
same translated source fibre
\(\mathcal F(u+s+\delta,\cdot)\).  It does not compare that fibre with
\(\mathcal F(u+s,\cdot)\).  Hence no unrecorded physical-coordinate
translation derivative or source-function regularity is being assumed.
Controlling that separate difference, and globally allocating it, remains
outside the theorem.

## 6. Compact-parameter extension: PASS

Let \(c\in I\Subset(0,\infty)\) and \(s/b\in J\Subset(0,\infty)\).
The tilted rate is determined uniquely by
\[
 {s\over b}={m(\gamma c)\over\gamma}
 =\frac1\gamma+\frac c{1+e^{\gamma c}}.
\]
The right side is strictly decreasing in \(\gamma>0\), and
\(1\le m(d)\le\mu\) gives
\[
 {1\over\max J}\le\gamma\le{\mu\over\min J}.
\]
Thus \(\theta=1-\gamma\) may be positive or negative, but remains in a
compact subset of \((-\infty,1)\).  The tilted exponential rates,
Bernoulli denominators, centered moments, characteristic-function decay,
and local-CLT normalizer are consequently uniform over \(I\times J\).
In particular the all-plus conditional density at its tilted mean stays
comparable to \(b^{-1/2}\).

Mean stationarity at \(\mu\) is not used here.  The \(c\)-score's linear
drift is removed by the same \(\tau\)-integration by parts, so its
coefficient may be \(O_{I,J}(1)\).  The raw \(s\)-derivative proof is
unchanged.  At each evaluation point the chosen tilt is held fixed;
integrating the uniform derivative bounds over a \(c\)- or \(s\)-path
introduces no \(\theta'(c)\) or \(\theta'(s)\) term.  Therefore the
\(C^1_c\) and \(C^1_s/b\) multiplier bounds, including joint labels,
extend uniformly to the compact family.

## Scope retained

The audit does not infer a common global source profile, a physical
modular comparison across receiver cones, or a Carleson bound for their
reuse.  Those are precisely the remaining steps before any general
strong-maximal conclusion.
