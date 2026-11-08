# Qualified angular budget bridge

## Result of this audit

**Important limitation.**  The angular mean \(B(w)\) need not enter the
active region \(B\ge1\) on any useful part of a high-efficiency physical
source.  Finite-\(N\) dilation can concentrate a bounded source into an
\(O(D^{-1})\) angular packet while preserving its endpoint ratio exactly;
then \(B_D(w)\to0\) uniformly.  See the dilation obstruction note in this
round.  Thus the estimates below are conditional local tools once an
independent active-region selection has been proved; they do not establish
such a selection from global budget or efficiency.

The source budget does rule out the previously feared generic
\(L^2\)-defect \(O(\sqrt\rho)\) on the active \(B\ge1\) radial region.
The correct supporting slope is the **right** slope of \(\Phi_N\), not
the left slope used only for one-sided gain events.  It gives a factor
\(N\) of angular quadratic rigidity:
\[
 {\mathbb E}_{\sigma,\eta}|F_\sigma(w,\eta)-B(w)|^2
 \le {C_A\,\mathcal D(w)\over N\Phi_N(B(w))}
 \qquad(B(w)\ge1).                                        \tag{1}
\]
Here \(B\) is the single common source radial profile obtained by angular
and label averaging, and \(\mathcal D\) is its actual angular Bregman
defect.

The conditional labelled first-response theorem and its finite-\(V\)
Dirichlet implementation are now proved locally.  On defect of natural
size they give an \(O(\rho/\sqrt N)\) action-comparison error, namely
\(O(\sqrt L/N)=o(L/N)\), after the exact radial-Jacobian cancellation
recorded below.  If an additional source-to-receiver localization confines
the relevant subunit component to a radial strip of probability
\(O(L^2/N)\), the bounded first-response estimate contributes only
\[
 O(\rho L^2/N)=o(L/N).                                    \tag{2}
\]

Even after an independent active-region selection, the global conclusion
is conditional on the still-unproved cone-packing/allocation statement:

1. a cone-packing/allocation bound which prevents the same source defect
   from being charged by many receiver cones.

The source budget alone implies neither statement.  In particular this
note does not claim a general strong-maximal upper bound.

## 1. One common source profile and its exact global budget

Here and below \(L=\log N\).  The unrelated finite physical depth is
denoted by \(D\); this avoids confusing the critical \(KL\) buffers with
the coarea depth.

The physical source radial domain is
\[
 \Omega_D=\{w\in\mathbb R^r:D\alpha_j+w_j>0,\ 1\le j\le r\}.
\]
This is the global domain.  It must not be replaced by the local
positive-orthant increment domain \(s\in\mathbb R_+^r\).

Let \(d_j=b_j-1\), \(v_j=D\alpha_j+w_j\), and let
\(F_\sigma(w,\eta)\in[0,A]\) be a raw terminal on orthant label
\(\sigma\) and product Dirichlet angle \(\eta\).  Average using uniform
probability over labels and angles:
\[
 B(w)={\mathbb E}_{\sigma,\eta}F_\sigma(w,\eta).           \tag{3}
\]
For a chosen supporting slope \(g_{B(w)}\), put
\[
 D_{B(w)}(x)=\Phi_N(x)-\Phi_N(B(w))-g_{B(w)}(x-B(w)),
\quad
 \mathcal D(w)={\mathbb E}_{\sigma,\eta}D_{B(w)}
 (F_\sigma(w,\eta)).                                      \tag{4}
\]
The affine term averages to zero at every source \(w\), so the exact
physical modular identity is
\[
 {\cal J}_D(F)-{\cal J}_D(B)
 =C_D\int_{\Omega_D} e^{-\sum w_j}
     \prod_jv_j^{d_j}\mathcal D(w)\,dw.                   \tag{5}
\]
This includes arbitrary orthant labels: conditional sign flips merely
permute the label sum.  It is an actual global source-budget identity,
not a family of local fibre budgets.

## 2. Right-slope Bregman rigidity on the active region

For \(N\ge2\) and \(B\in[1,A]\), use the right slope
\(g_B=\Phi_N'(B+)\).  There is a constant \(c_A>0\), independent of
\(N\), such that
\[
 \Phi_N(x)-\Phi_N(B)-\Phi_N'(B+)(x-B)
 \ge c_A N\Phi_N(B)(x-B)^2
 \quad(0\le x\le A).                                      \tag{6}
\]

To check (6), on \(h>1\)
\[
 \Phi_N'(h)=(1+\log h)^{N-2}(N+\log h),
\]
and
\[
 \Phi_N''(h)
 ={(N-1)(N-1+\log h)(1+\log h)^{N-3}\over h}.
\]
On a \(1/N\)-neighborhood of \(B\), Taylor's theorem gives a lower
bound of order \(N^2\Phi_N(B)(x-B)^2\), stronger than (6).  Outside that
neighborhood, monotonicity of \(\Phi_N'\) supplies a gap of order
\(N\Phi_N(B)|x-B|\), which implies (6) because \(x,B\in[0,A]\).  If
\(x\le1\), use \(\Phi_N(x)=x\) and the same tangent-gap bound at one.
The constants can be made uniform over the fixed cap interval.

Since the affine term in (4) has zero angular-label mean, (6) proves
(1).  This is the structural improvement over a bare \(L^2\) estimate:
a rare order-one low value can only be balanced by a positive fluctuation
whose endpoint modular cost is order \(N\) times its squared size.

For comparison, at the literal corner \(B=1\), the left slope gives the
one-sided inequality
\[
 D_1(x)=\Phi_N(x)-x\ge(N-1)(x-1)_+ .                       \tag{7}
\]
If the angular-label mean is one, then
\[
 {\mathbb E}|F-1|\le {2\mathcal D\over N-1},\qquad
 {\mathbb E}|F-1|^2\le {C_A\mathcal D\over N-1}.           \tag{8}
\]
Equation (8) is useful for deficits; (6) is the stronger symmetric
estimate needed for first-response control.

The exact corner expansion is
\[
 D_1(1+t/N)=\Phi_N(1+t/N)-1-t/N=e^t-1+o(1).                \tag{9}
\]
No \(e^t-1-t\) coefficient is used here.

## 3. Proved conditional finite-depth interface

For one block and one action coordinate, let \(K_c^{V,s}\) be the
conditional angular kernel with the complete increment and orthant-flip
law retained.  The conditional flux theorem gives the dimension-free
derivative on the tangent hyperplane.  At finite depth, zero-extending
\(G=F-B\) from the source simplex gives the exact Dirichlet-scaled form
\[
 \|(K_c^{V,s}-K_{c_0}^{V,s})G\|_{L^2(d\sigma_\xi)}
 \le C|c-c_0|
 \left({V+s\over V}\right)^{d/2}
 \|G\|_{L^2(d\sigma_\eta)}.                               \tag{10}
\]
The label norm on both sides is uniform Haar measure on the orthant label
group, and the angle norm is uniform Dirichlet measure; the conditional
jump/sign law belongs to the operator, not to the source norm.  This is
not an unweighted simplex inequality.  Multiplication by
\(V^d\) in squared norm changes the factor on the right exactly to
\((V+s)^d\), the source radial Jacobian.
The constant \(B(w)\) itself maps to \(B(w)\) under every action; only
\(F-B\), not \(B\), is zero-extended.  For fixed \(r\), the same statement
holds for a vector action varied in one coordinate, or for a common action
by summing its fixed number of block derivatives.  In the product
notation of (11), the independent block Jacobian factors multiply.  It
makes no claim about a multi-parameter maximal operator.

On \(B\ge1\), insert \(G=F-B\), (1), and
\(|c-c_0|=O(\rho)\).  This yields
\[
 \left(\prod_jV_j^{d_j}\right)^{1/2}
 \|(K_c^{V,s}-K_{c_0}^{V,s})(F-B)\|_2
 \le C\rho
 \left(\prod_j(V_j+s_j)^{d_j}\right)^{1/2}
 \sqrt{\frac{\mathcal D(w)}{N\Phi_N(B(w))}}.              \tag{11}
\]
Thus the local response scale is
\(\rho/\sqrt N=\sqrt L/N=o(L/N)\), before the separate
receiver-cone allocation.

## 4. The lower core and transition strip

Take a normalized critical source budget with radial threshold
\(q_0=N\mu-d_N\).  Since \(\Phi_N(F)\ge F\), on the far lower source core
\(\sum_jw_j\le q_0-KL\) the actual global modular gives
\[
 \int_{\rm low}F\,d\lambda_{\rm src}
 \le e^{q_0-KL}{\cal J}_D(F)=O(N^{-K})                    \tag{12}
\]
after the usual normalization.  This is a source statement and remains
valid with labels and angles summed, where
\[
 d\lambda_{\rm src}
 =C_D\prod_jv_j^{d_j}\,dw\,d\sigma(\eta)\,d\mathbb P(\sigma)
\]
is the **unweighted** coarea measure.  The physical Lebesgue modular has
the additional factor \(e^{-\sum w_j}\).  It removes the potentially
uncontrolled \(B<1\) angular variance from the far core.

To continue from (12), one needs the intended critical-envelope
localization: apart from the suppressed lower core, the subunit part which
can affect a near-success receiver is confined to a radial strip of width
\(O(L)\).  This confinement is an additional source-to-receiver statement,
not a consequence of (12) by itself.  Once it is available, its critical
block-law probability is
\[
 O\bigl(L\,f_{N,\mu}(q_0+O(L))\bigr)=O(L^2/N).              \tag{13}
\]
Provided the conditional first-response theorem supplies a bounded-test
first derivative, the contribution of this strip over an action window
of length \(O(\rho)\) is at most
\[
 C\rho\,O(L^2/N)=o(L/N),                                  \tag{14}
\]
because the terminal is capped.  This avoids treating a whole lower tail
as a probability-\(O(\rho)\) set with arbitrary \(L^2\) defect.

Equations (12)--(14) still require transport through the actual
source-to-receiver coupling.  They are not a licence to apply a radial
theorem fibre by fibre.

## 5. Exact cone allocation that is still missing

For a receiver radial set \({\cal U}\), Fubini applied to the nonnegative
defect density in (5) gives the exact cone identity
\[
 C_D\int_{\cal U}\int_{s\ge0}
 e^{-\sum(u_j+s_j)}
 \prod_j(D\alpha_j+u_j+s_j)^{d_j}\mathcal D(u+s)\,ds\,du
 =\int_{\Omega_D}M_{\cal U}(w)\,d\Lambda_{\rm def}(w),     \tag{15}
\]
where
\[
 d\Lambda_{\rm def}
 =C_D e^{-\sum w_j}\prod_jv_j^{d_j}\mathcal D(w)\,dw,
\qquad
 M_{\cal U}(w)=|\{u\in{\cal U}:-D\alpha<u\le w\}|.         \tag{16}
\]
This is the same extra upper-cone multiplicity as in the global
calculation of angular_budget64.  A global budget controls
\(\int d\Lambda_{\rm def}\), but cannot control
\(\int M_{\cal U}d\Lambda_{\rm def}\): defect may be concentrated where
\(M_{\cal U}\) is large.

The following is only a prospective sufficient physical packing statement,
not a consequence of the defect identity:
\[
 \int M_{\cal U}\,d\Lambda_{\rm def}
 \le C_{\rm pack}\,\nu({\cal U}),                          \tag{17}
\]
with the appropriate policy likelihood factor included when a
receiver-dependent action is selected.  If such a weighted constant were
\(O(1)\), (7) gives average deficit \(O(1/N)\), and (11) gives average
first-response error \(O(\sqrt L/N)\); both are \(o(L/N)\).

Neither (17) nor an acceptable replacement is implied by modular budget,
Jensen, or the fixed-s Dirichlet cancellation.  In particular no automatic
defect-per-receiver bound is being inferred.  This is the exact remaining
common-source compatibility and receiver-cone allocation gap.

## 6. Non-counterexample check

For the two-level angular spike of mass \(m\), height \(1+t/N\), and
mean-one compensation below one, the defect is
\[
 m(e^t-1+o(1)).
\]
Taking \(m=1/\widetilde Z_s\) is compatible with a policy-envelope
estimate and has only a fixed-\(t\) \(O(1/N)\) received perturbation.
It is repairable in the fixed radial theorem.  It therefore does not
disprove (10), (17), or a qualified \(o(L/N)\) transfer.
