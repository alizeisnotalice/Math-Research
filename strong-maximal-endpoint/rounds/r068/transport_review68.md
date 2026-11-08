# Band-averaged transport audit

## Verdict

The proposed change of variables eliminates the **graph trace** for a
nonnegative band-averaged source transport.  It is rigorous provided the
conditional kernels satisfy the stated strong \(L^1\) envelope in the
radial parameter.  It does not by itself control the original pointwise
fixed-\(s\) graph operator, produce a reference-action success statement,
or prove a weak bound by \(Q_N\).

Use global logarithmic centroid coordinates
\[
 w=\sum_{i=1}^bv_i,\qquad X=v-\frac wb{\bf1}\in H_b.
\]
With Hausdorff Lebesgue measure \(dH(X)\) on \(H_b\),
\[
 dv=\frac1{\sqrt b}\,dw\,dH(X),\qquad
 dx=\frac{2^b}{\sqrt b}e^{-w}\,dw\,dH(X)
\]
after summing physical orthants.  Let \(d\lambda\) be normalized Haar
measure on the finite orthant-label group.  For \(G\ge0\), define
\[
 A_u(X,\sigma)=\int_I e^{-(u+s)}
       [K_c^sG_{u+s}](X,\sigma)\,ds,
 \qquad G_w(X,\sigma)=G(w,X,\sigma).                     \tag{1}
\]
The physical source modular in these coordinates is
\[
 {\cal J}(G)=\frac{2^b}{\sqrt b}
 \int_{\mathbb R}e^{-w}\|G_w\|_{L^1(dH\,d\lambda)}\,dw.  \tag{2}
\]
There is no \(V^{b-1}\) factor in this global coarea formula.  Such a
factor belongs only to a normalized simplex-angle chart, not to the
translation-invariant \(H_b\) convolution coordinates used here.

## 1. Exact band transport calculation

With \(w=u+s\),
\[
 A_u(X,\sigma)=\int_{\mathbb R}{\bf1}_{\{w-u\in I\}}\,
 e^{-w}[K_c^{\,w-u}G_w](X,\sigma)\,dw.
\]
Positivity gives the pointwise bound
\[
 \sup_u A_u(X,\sigma)
 \le\int_{\mathbb R}e^{-w}
       \sup_{s\in I}[K_c^sG_w](X,\sigma)\,dw.             \tag{3}
\]
There is no multiplicity factor in (3): for each fixed \(u,w\), exactly
one value \(s=w-u\) occurs, and the indicator is at most one.

The needed radial conditional-kernel input is the exact envelope identity
of envelope_proof68:
\[
 \left\|\sup_{s\in I}K_c^s g\right\|_1
 \le Z_I\|g\|_1,\qquad
 Z_I=1+\log\frac{H_{b,c}(s_1)}{H_{b,c}(s_0)},\quad
 H_{b,c}(s)=\sum_{k=0}^b{b\choose k}
 { (s-ck)_+^{b-1}\over(b-1)!}.                            \tag{4}
\]
Then Tonelli applied to (3) gives the advertised trace-free bound
\[
 \frac{2^b}{\sqrt b}\int\sup_u A_u(X,\sigma)\,dH(X)\,d\lambda(\sigma)
 \le Z_I\,{\cal J}(G).                                   \tag{5}
\]
This applies equally to a measurable choice \(u=u(X,\sigma)\), since
\(A_{u(X,\sigma)}(X,\sigma)\le\sup_uA_u(X,\sigma)\).

The proof of (4) is not a total-variation differentiation argument.
The conditional kernel supports grow by one-sided nested shells as \(s\)
increases.  The new shell mass is
\(\partial_s\log H_{b,c}(s)\,ds\), so integrating only the new portions
gives
\[
 \int\sup_{s\in I}k_{c,s}(y,J)\,dy\,d\#J
 =1+\log{H_{b,c}(s_1)\over H_{b,c}(s_0)}.
\]
In contrast, the full distributional derivative has total-variation
mass \(2\,\partial_s\log H_{b,c}(s)\); using
\(\int\|\partial_sK_c^s\|_{1\to1}ds\) would give the wrong factor and is
not an acceptable proof of (4).

The measure convention is essential but harmless: \(dy\) is Lebesgue
measure in the convolution coordinate on \(H_b\), and \(d\#J\) is
counting measure on labels.  Replacing these by \(dH\,d\lambda\)
rescales the kernel and its reference measure together, so the normalized
envelope mass is still \(Z_I\); the sole physical coarea factor is then
the \(2^b/\sqrt b\) displayed in (5).  For \(b=2\),
\(H_{2,c}\) is positive and piecewise linear on \((0,\infty)\), hence
absolutely continuous; \(\log H_{2,c}\) is absolutely continuous on every
compact \(I\Subset(0,\infty)\), including bands whose endpoints or
interior meet a knot \(ck\).  Thus no endpoint atom or extra constant is
missing in (4).

## 2. What has been gained

For a fixed \(s\), the graph operator has the nonconvolution form
\[
 X\mapsto K_c^s\!\left[G_{u_*(X)+s}\right](X),
\]
whose input changes with the output \(X\).  A fixed-function \(L^2\)
bound cannot be applied to it.

In the \(w=u+s\) representation, each source fibre \(G_w\) is
first acted on by a genuine fixed-function convolution \(K_c^s\), and
only then is the \(s\)-supremum taken.  Formula (5) is therefore a real
strong \(L^1\) replacement for the missing graph trace.  It is compatible
with arbitrary measurable \(u(X)\).

## 3. What it does not give for an actual halo

The physical receiving average has radial law
\[
 {\cal R}_cG(u,X)
 =\int f_{S,c}(s)\,K_c^sG_{u+s}(X)\,ds,                   \tag{7}
\]
with the appropriate labelled version.  The weight in (1) is
\(e^{-(u+s)}=e^{-w}\), a source modular weight, not automatically
\(f_{S,c}(s)\).  Writing
\[
 f_{S,c}(s)=e^{-s}H_c(s)
\]
shows the missing comparison explicitly:
\[
 e^{-u}{\cal R}_cG(u,X)
 =\int e^{-(u+s)}H_c(s)K_c^sG_{u+s}(X)\,ds.               \tag{8}
\]
Hence a pointwise or band-integrated success certificate for (7) yields
a lower bound for \(A_u\) only after controlling the band mass and the
size of \(H_c\) on that band.  Success concentrated at one value of \(s\)
does not lower-bound the integral in (1).

There is a useful nonmonotone thickness reduction.  For each angular-label
state, let \(u_*(X)\) be the infimum of its successful radial set.  Without
assuming that the halo is monotone in \(u\),
\[
 \int_{\{u:\ {\rm successful}\}}e^{-u}\,du\le e^{-u_*(X)},
\]
because every successful \(u\) is at least its infimum.  Thus one selected
success point can in principle control the whole positive radial halfline.
To use this rigorously one still needs a finite lower radial cutoff and a
measurable selection/approximation of \(u_*(X)\) (for example via rational
levels and a limiting argument).  This is a selection issue, not a reason
to require halo monotonicity.

## 4. No direct fixed-reference conclusion

The estimate (5) is positive and \(L^1\); it compares neither actions nor
the full radial averages.  It does not imply
\[
 {\cal R}_cG>1\quad\Longrightarrow\quad
 {\cal R}_{c_0}G>1-\varepsilon
\]
for a fixed \(c_0\), nor a \(Q_N\) weak estimate.  Such a conclusion
would additionally require:

1. a band-capture lemma transferring a definite fraction of (7) into a
   lower bound for (1), with the \(H_c\) factors recorded;
2. an action-comparison estimate on the same band for a fixed source
   fibre, using the conditional \(C^1\) theorem and a source defect
   budget; and
3. a radial stopping/packing lemma converting the selected graph or
   intervals into physical halo volume.

The new transport is therefore a useful candidate for item 3: it removes
the fixed-\(s\) trace obstruction after band averaging.  It does not
remove the independent source-budget, action-sharing, or cone-packing
requirements, and supplies no general-\(C_N\) counterexample.

Finally, the bridge has the expected exact radial scaling.  Under
\(G^D(w,X,\sigma)=G(w-bD,X,\sigma)\), a direct change of variable gives
\(A_u^D=e^{-bD}A_{u-bD}\), while the physical modular in (2) is multiplied
by \(e^{-bD}\).  Hence (5), including its coarea constant, is scale
covariant.  This bookkeeping does not create the still-missing comparison
from a variable-action success point to a fixed-reference \(A_{u,c_0}\).

## Appendix: audit of the fixed-block mixed obstruction

Final independent audit: PASS after the direct-BV repair in the final mixed proof.

The conclusion of `mixed_obstruction68.md` is valid for fixed
\(b\ge4\), \(c>0\), and \(s>c\) away from the nodes \(kc\).  Its
displayed (7)--(10) now supply the necessary direct one-dimensional BV
argument, rather than resting on a radial first-response theorem
formulated for another Fourier class.

For \(\xi(q)=q(b-1,-1,\ldots,-1)\), the face which fixes the first
coordinate at \(c\) is exactly
\(E_q=e^{iq(bc-s)}\).  Every other face, and the full multiplier, has a
one-dimensional representation with oscillation frequency \(bq\).  On a
small node-free \(s\)-interval their zero-extended marginal densities
\(g_s\), as well as \(\partial_sg_s\), are BV (in particular the latter
are uniformly \(L^1\)).  For the full multiplier this follows from
\[
 g_s(z)={q_c(z)f_{b-1,c}(s-z)\over f_{b,c}(s)},
         {\bf1}_{[0,s]}(z).
\]
For a low face use the analogous formula with \(b-1\) conditioned
coordinates.  The endpoint terms cause no delta: when \(b\ge4\), the
relevant spline vanishes to sufficient order at its moving endpoint;
the finitely many interior spline knots only give BV jumps.  This is why
the stated fixed-\(b\ge4\) range is adequate.

Integration by parts gives \(\Phi=O(q^{-1})\), while differentiating its
one-dimensional Fourier integral gives
\(\partial_s\Phi=-iq\Phi+e^{-iqs}\widehat{\partial_sg_s}(bq)=O(1)\), and the
same bound holds for each low face.  Therefore, writing the mean face as
\(b^{-1}E_q+r_q\), one has \(r_q=O(q^{-1})\) and
\(\partial_sr_q=O(1)\).  Differentiating the exact flux identity gives
\[
 \partial_s\partial_c\Phi
 =iq\,{A_{b,c}(s)\over b}E_q+O(1).
\]
The sign is positive, and \(A_{b,c}(s)>0\).  Hence the unboundedness
claim follows.  A differentiation of the bare \(O(q^{-1})\) remainder
in the face expansion would not be legitimate; the separate BV derivative
bound is the essential point.

## Appendix: audit of the radial envelope and rate

`radial_envelope68.md` has the correct normalization: with Hausdorff
measure on \(H_b\), the conditional density is
\((\sqrt b H(s))^{-1}{\bf1}_{D_{s,J}}\), and summing labels gives
\(\sum_J|D_{s,J}|=\sqrt bH(s)\).  Thus nested supports give exactly
\(1+\log(H(d)/H(a))\), including across \(b=2\) knots; total variation
would instead be \(2\partial_s\log H\).  The approximate-identity
argument establishes sharpness for this positive maximal convolution,
not a halo sharpness statement.

The compact-parameter rate \(B_b(c,s)=\gamma+O(b^{-1})\) is also
consistent: the tilt \(\theta=1-\gamma\) has mean
\(\gamma^{-1}+c/(1+e^{\gamma c})\), the tilt family has uniformly
nondegenerate variance and a uniform nonlattice Fourier bound, and the
derivative of the tilted density at its mean is \(O(b^{-3/2})\) while the
density is \(\asymp b^{-1/2}\).  This yields the required logarithmic
derivative error \(O(b^{-1})\).  Its stated scope must remain compact
positive ranges of \(c\) and \(s/b\); it is not a uniform assertion down
to \(s=0\), across unbounded parameters, or for a growing product-block
count.
