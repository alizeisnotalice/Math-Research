# Mathematical assertion and scope

Let D = {(u,v): |u|+|v| <= 1}. For 0 < alpha < 2 use the coordinate jump generator

\[
Af(x)=\sum_{j=1}^2\int_{\mathbb R}
[f(x+t e_j)+f(x-t e_j)-2f(x)]|t|^{-1-\alpha}\,dt.
\]

This is a positive multiple of the operator with Fourier symbol `-(|xi1|^alpha+|xi2|^alpha)`, so the sign condition is unchanged by that normalization.

The Lin-Tomar comparison principle requires an even, integrable K with K >= 0 a.e., K >= 1 on D a.e., and AK >= 0 **as a distribution outside the origin**. It gives `C2 <= (1/2) integral K`. The invertible map `(u,v) -> (u+v,u-v)` identifies the diamond and centred square constants; its Jacobian multiplies input mass and receiving area by the same factor. Appendix A of the paper extends the transfer from smooth compact inputs to bounded compact inputs and then all L1 inputs by mollification and monotone truncation, preserving all positive scales.

The original and new certificates both have alpha=47/40, R=33/20 and h=1/32. They consist of a truncated singular power plus finitely many cubic tensor splines. The singularity r^(-alpha) is integrable in dimension two because alpha<2. Lemma 3.1 turns the original a.e. pointwise C4 into the required distributional statement away from zero. No integrability of AK near the origin, and no representation as a finite signed measure minus an origin atom, is assumed.

The original complete checker run certifies C1-C4. Outside the entire kernel support K=0, C1 implies AK>=0 a.e. by the nonnegative jump integrals. The spline can be nonzero beyond r=R: the old C4 report includes that region. Axes and r=R form null sets; the paper proves local integrability of the difference integral away from the origin, which is the required bridge to distributions.

Set `kappa=10000993/10000000` and Knew=K/kappa. The verifier checks exact rational equality for c and every coefficient class, with all other mathematical fields unchanged. It actually proves Knew>=1 on D with the original checker at 192-bit precision; it also independently rechecks the stronger K>=kappa statement. These comparisons use `lower().fmpq()` exact endpoints, not decimal displays.

Positive scaling preserves evenness, integrability and nonnegativity, and for every nonnegative smooth test phi supported away from zero,

\[
\int K_{new}A\phi=\kappa^{-1}\int KA\phi\ge0.
\]

Thus new C4 follows from the recorded original full C4 and linearity. The default verifier does not recompute the old leaf bounds. Fresh full-check commands are supplied in the README for readers who want to replace that historical-run premise by a new numerical execution.

Direct integration gives

\[
\frac12\int K_{new}
=\frac{c_{new}\alpha R^{2-\alpha}}{2-\alpha}
 +\frac{h^2}{2}\sum_{i,j}a^{new}_{ij}.
\]

The sum is independently expanded over all indices. Arb evaluates the remaining positive real power with outward enclosure; exact rational endpoints prove

\[
3.590867754338446537
<\frac12\int K_{new}
<3.590867754338446538
<\frac{3590867754339}{10^{12}}.
\]

The comparison principle therefore gives the strict bound `C2 < 3.590867754339` for general L1 inputs and all scales. The strict and non-strict level conventions yield the same global weak constant by a threshold limit; no pointwise equality of the two level sets is claimed.

The new coefficients are general rationals, not necessarily dyadic. The public certificate schema and the unchanged checker both accept rational strings. No claim is made that the old Lean data encoding accepts the new coefficients without further work, or that this package performs a Lean build.

This proves a valid upper bound, not optimality of the kernel family, equality with C2, or existence/nonexistence of an extremizer.
