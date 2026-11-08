# Finite-depth dilation obstruction to the global angular-mean baseline

## Statement

The global angular-label average
\[
 B(w)={\mathbb E}_{\sigma,\eta}F_\sigma(w,\eta)
\]
cannot be forced into the active range \(B\ge1\) by endpoint efficiency,
bounded height, or the physical modular budget.  The obstruction is exact
at every fixed finite \(N\): ordinary isotropic dilation preserves the
strong-maximal endpoint ratio while shrinking the angular support of a
source packet to the barycenter at rate \(D^{-(N-r)}\).

This is not a counterexample to a strong-maximal inequality.  In
particular, an everywhere subunit input has no receiving halo.  It
obstructs only the proposed mechanism that tries to obtain a global
action-independent radial baseline by averaging over source angles and
then invokes the \(B\ge1\) Bregman rigidity everywhere.

## 1. Exact dilation invariance

For a bounded measurable \(f\) on \(\mathbb R^N\), define
\[
 f_D(x)=f(e^D x),\qquad D>0.                              \tag{1}
\]
The strong maximal operator commutes with this dilation:
\[
 M_Sf_D(x)=M_Sf(e^Dx).                                    \tag{2}
\]
Indeed rectangles are carried bijectively to rectangles by
\(R\mapsto e^D R\), with the normalized averages unchanged.  Therefore,
for every threshold \(\lambda\),
\[
 |\{M_Sf_D>\lambda\}|=e^{-ND}|\{M_Sf>\lambda\}|.          \tag{3}
\]
The endpoint modular has the same scaling,
\[
 \int_{\mathbb R^N}\Phi_N(f_D(x))\,dx
 =e^{-ND}\int_{\mathbb R^N}\Phi_N(f(y))\,dy.               \tag{4}
\]
Thus every endpoint volume-to-modular ratio is exactly invariant.  The
cap is unchanged.  A strict rectangle witness for \(f\) maps to a strict
rectangle witness for \(f_D\), with the same strict margin.

## 2. Angular concentration at finite depth

Assume the source is supported away from coordinate planes in a compact
logarithmic box.  More precisely, in a fixed orthant put
\[
 v_i=-\log|y_i|,\qquad a_i\le v_i\le b_i \quad(1\le i\le N), \tag{5}
\]
with fixed finite \(a_i,b_i\).  For \(x=e^{-D}y\),
\[
 -\log|x_i|=D+v_i.                                       \tag{6}
\]
Partition the coordinates into fixed blocks of sizes \(b_j\), write
\[
 S_j=\sum_{i\in{\cal B}_j}v_i,\qquad
 V_j=\sum_{i\in{\cal B}_j}(-\log|x_i|)
 =b_jD+S_j.                                               \tag{7}
\]
The within-block source angle is
\[
 \eta_i={D+v_i\over b_jD+S_j}\qquad(i\in{\cal B}_j).
\]
It obeys the exact identity
\[
 \eta_i-\frac1{b_j}
 =\frac{b_jv_i-S_j}{b_j(b_jD+S_j)}.                       \tag{8}
\]
The numerator is uniformly bounded on the box (5), while the denominator
is bounded below by a positive multiple of \(D\).  Hence the angular
support in each \((b_j-1)\)-simplex lies in an \(O(D^{-1})\) neighborhood
of its barycenter.  Its uniform Dirichlet probability is at most
\[
 C_{N,a,b}D^{-\sum_j(b_j-1)}
 =C_{N,a,b}D^{-(N-r)}.                                   \tag{9}
\]
This bound is uniform in the radial source coordinate: outside the
radial projection of the support it is zero, and inside it (8) applies.

Let \(F_{D,\sigma}(w,\eta)\) denote the exact finite-depth coarea terminal
of \(f_D\), including its orthant labels.  From \(0\le F_D\le A\), (9),
and averaging uniformly over labels and angles,
\[
 0\le B_D(w)
 :={\mathbb E}_{\sigma,\eta}F_{D,\sigma}(w,\eta)
 \le C_{N,a,b}A D^{-(N-r)}                               \tag{10}
\]
for every \(w\) in the finite-depth global radial domain.  In particular
\[
 \sup_wB_D(w)\longrightarrow0\qquad(D\to\infty).          \tag{11}
\]
Equations (6)--(11) concern the physical finite-depth coarea variables.
They are not a statement about a deep-limit radial terminal defined only
on a local positive orthant.

## 3. Log-box exhaustion and cap-specific near-endpoint applicability

For every nonnegative \(f\) with finite \(\Phi_N\)-modular, define
\[
 f_j(x)=\min(f(x),j)\,
 {\bf1}_{\{\,1/j\le |x_i|\le j\ \text{for every }i\,\}}(x). \tag{12}
\]
Then \(f_j\uparrow f\) almost everywhere.  For every rectangle its average
increases to the average of \(f\), and the two suprema commute:
\[
 M_Sf_j\uparrow M_Sf.
\]
Consequently
\[
 |\{M_Sf_j>1\}|\uparrow|\{M_Sf>1\}|,\qquad
 \int\Phi_N(f_j)\uparrow\int\Phi_N(f).                    \tag{13}
\]
Whenever the limiting modular is positive, the corresponding endpoint
ratios converge.  Thus compact logarithmic-box inputs are sufficient for
the full finite-\(N\) endpoint problem, without a compact witness-family
or uniform-margin hypothesis.

For a fixed cap \(A\), start in the cap class \(0\le f\le A\).  The
additional amplitude truncation in (12) is then eventually irrelevant,
and the spatial cuts preserve the cap.  Hence the cap-specific endpoint
supremum can be approximated by capped compact-log-box inputs.  Dilation
of a chosen approximant preserves its ratio exactly by (3)--(4), and a
subsequent choice of large \(D\) makes its canonical global angular mean
strictly below one everywhere by (10).

In particular, for fixed \(N\) and fixed \(A\), a near-\(Q_N^{[A]}\)
input may be chosen separately from this compact class and then dilated to
produce a comparably efficient capped packet with \(B_D<1\) everywhere.
This is not a statement uniform in \(A\) or in \(N\), and it does not
identify the unrestricted \(Q_N\) as \(N\) varies.  The explicit
Jensen-sharp box in Section 5 is a different example: it is not asserted
to be near-\(Q_N^{[A]}\).

The obstruction is not removed by a literal outer-support normalization.
Choose \(D\) so large that the dilated compact main packet \(f_D\) is
supported in \(\{|x_i|<1/2\ \forall i\}\) and has
\(B_D\le1/4\).  Let
\[
 E=\{x:3/4\le|x_i|\le1\ \text{for every }i\},\qquad
 |E|=2^{-N},
\]
and, for \(\eta>0\), set
\[
 \varepsilon=\min\left\{\frac14,\,
       {\eta\,{\cal J}(f_D)\over |E|}\right\},\qquad
 \widetilde f_D=f_D+\varepsilon{\bf1}_E.                  \tag{14}
\]
The supports are disjoint.  Since \(\Phi_N(\varepsilon)=\varepsilon\),
\[
 {\cal J}(\widetilde f_D)\le(1+\eta){\cal J}(f_D),\qquad
 |\{M_S\widetilde f_D>1\}|\ge|\{M_Sf_D>1\}|.              \tag{15}
\]
The cap is preserved for every fixed \(A>1\), every coordinate has
essential outer radius exactly one, and the canonical global angular mean
obeys
\[
 \widetilde B_D\le B_D+\varepsilon\le\frac12.              \tag{16}
\]
Thus a cap-specific near-endpoint ratio loses at most the factor
\((1+\eta)^{-1}\), with \(\eta\) arbitrary, even after this elementary
bounding-box normalization.

This concerns a fixed-reference angular chart and its plain support
normalization only.  It does not exclude an energy-based effective-support
selection, a localized chart, or an adaptive chart; any such mechanism
would need a separate proof.

## 4. Why global angular Bregman rigidity cannot start

Take any nonconstant angular packet \(G_D(\eta)\in[0,a]\) with \(a<1\),
supported in the set from (9), and define a source terminal independent
of \(w\) on a compact radial packet.  For all sufficiently large \(D\),
its angular mean \(B_D\) is below one.  On the entire subunit interval,
\(\Phi_N(h)=h\), so with the only slope \(g=1\),
\[
 D_{B_D}(G_D)=\Phi_N(G_D)-\Phi_N(B_D)-(G_D-B_D)=0.          \tag{17}
\]
Nevertheless \(G_D\) has nonzero angular variance.  Therefore neither
the modular budget nor the defect can control angular \(L^2\) variation
below one.

For amplitudes exceeding one, the same dilation still makes
\(B_D<1\) everywhere, while the source has nonzero modular defect on its
tiny angular packet.  By (3)--(4), it can retain whatever bounded
near-endpoint ratio the undilated strict packet had.  Thus near-endpoint
efficiency does not repair the failure of global active-region coverage.

The right-slope inequality for a preselected local region with \(B\ge1\)
remains valid and useful.  What fails is the additional claim that the
canonical global angular mean supplies such a region at the needed
endpoint mass.

## 5. Explicit capped cube packet

There is a concrete finite-\(N\) packet illustrating the obstruction
without invoking any near-\(Q_N\) approximation.  Let
\[
 f_N(y)=a_N{\bf1}_{\{e^{-1}\le|y_i|\le1\ (1\le i\le N)\}}(y),\qquad
 a_N=1+N^{-2},\qquad N>r,                                 \tag{18}
\]
and dilate it as in (1).  It is capped, has the same terminal on every
orthant label, and each original orthant cube is a strict receiving witness
because \(f_N=a_N>1\) there.

For a block of size \(b\), let \(S\in[0,b]\) be the sum of the original
logarithmic coordinates \(v_i\in[0,1]\), and let
\(f_{{\rm IH},b}\) be the Irwin--Hall density.  The exact uniform
Dirichlet angular average of the dilated block indicator is
\[
 {\Gamma(b)f_{{\rm IH},b}(S)\over(bD+S)^{b-1}}.
\]
Because every orthant terminal is identical, the uniform all-label angular
mean of the full \(r\)-block packet is exactly
\[
 B_D(w)=a_N\prod_{j=1}^r
 {\Gamma(b_j)f_{{\rm IH},b_j}(S_j)\over(b_jD+S_j)^{b_j-1}}.
                                                                    \tag{19}
\]
Since \(f_{{\rm IH},b}\le1\) and
\(\Gamma(b)\le b^{b-1}\),
\[
 B_D(w)\le a_ND^{-(N-r)}.                                 \tag{20}
\]
Taking \(D=2\), the right side is below one for every \(N>r\).  Hence the
angularly radialized input \(B_D\) has no receiving halo, although the
original physical input has strict cube witnesses in every orthant.

The comparison is also sharp at the modular level.  Angular averaging
preserves the ordinary integral (the common \(2^N\) orthant factor cancels),
while \(B_D<1\) makes
\(\Phi_N(B_D)=B_D\).  Thus
\[
 {{\cal J}(f_{N,D})\over{\cal J}(B_D)}
 ={\,\Phi_N(a_N)\over a_N}
 =(1+\log a_N)^{N-1}\longrightarrow1.                    \tag{21}
\]
This disproves any sharp Jensen-budget-plus-near-unit-receiving argument
that would force the canonical angular mean into the active region.
It is not a near-\(Q_N\) or strong-maximal counterexample: its role is
only to isolate the failure of that baseline mechanism.  The cap-specific
log-box exhaustion in Section 3 is the separate near-endpoint route.
