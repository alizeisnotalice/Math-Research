# Independent analytic review of108

Verdict: the stated peak theorem and conditional104 coefficient transfer
are valid, subject to the explicitly accepted098/104 dependencies below.
This is not a global maximal-operator bound. The proof here is a review
of the actual pinned theorem, not an appeal to a general BL result whose
full source has not been read.

## Heat flow and all limits

Take unit u_i and positive c_i with sum c_i u_i u_i^T=I_r. Trace gives
sum c_i=r. The only required functions are exp(-alpha_i|s|). Their
truncations on[-R,R] are bounded, nonzero and integrable. Heat convolution
is strictly positive for every t>0, so logarithmic derivatives exist even
when a truncated function has jumps at the endpoints.

On a<=t<=b, convolution is bounded by C exp(-s^2/(8b)). Derivative ratios
are expectations of (y-s)/(2t) and (y-s)^2/(4t^2)-1/(2t) under the positive
heat-source measure supported on[-R,R]. Thus the stated linear/quadratic
bounds are valid. In the product h, the frame identity converts sum
c_i(u_i.x)^2 into |x|^2. Time derivatives, gradients and the divergence
used in the proof are a polynomial times C exp(-|x|^2/(8b)). These bounds
justify differentiation under the integral, and the smooth radial cutoff
has vanishing boundary integral as its radius tends to infinity.

With ell_i=log g_i,

 F'(t)=int h[sum c_i ell_i'^2-|sum c_i u_i ell_i'|^2]>=0.

The inequality follows from AA^T=I for A's columns sqrt(c_i)u_i;
A is a contraction on its entire coefficient space. It does not assume
that the projected directions are mutually orthogonal.

For t tending to infinity, x=sqrt(t)y and sum c_i=r cancel every power
of t. The inequality (s-z/sqrt(t))^2>=s^2/2-z^2/t, with |z|<=R and t>=1,
gives the integrable majorant C exp(-|y|^2/8). The limiting Gaussian
integral is exactly product (int f_i,R)^c_i. For t decreasing to0,
pointwise heat approximation holds away from finitely many hyperplanes
u_i.x=+-R; these have Lebesgue measure0 because every u_i is unit.
Fatou has the needed direction: initial integral<=liminf F(t). This
initial integral is bounded by the large-time limit by monotonicity.
Finally R tending to infinity increases each side, so monotone convergence
proves the needed geometric BL inequality with constant1. No finite
measure approximation is left unjustified, and equality cases are unused.

## Original projected law, zeros and moments

The fixed symbol is psi(xi)=sum_{v=(1,+-1,...,+-1)}|v.xi|, with32 terms,
unit Cauchy weights, Fourier inverse (2pi)^(-r), orthonormal quotient
Lebesgue measure and common time. Sign parity gives sum vv^T=32I_6.
For every orthogonal rank-r projection P, a_v=Pv and z_v=|a_v|^2 obey
sum a_v a_v^T=32I_E and sum z_v=32r. In particular positive projected
columns span E. For z_v>0 the choices u_v=a_v/sqrt(z_v), c_v=z_v/32 and
f_v(s)=exp(-32|s|/sqrt(z_v)) satisfy the BL hypotheses and reproduce the
symbol exactly. Its one-dimensional mass is sqrt(z_v)/16. When z_v=0
that original symbol term vanishes identically and its index is omitted;
no undefined direction, positive weight or0^0 is used. The convention
z log z=0 is only its continuous extension in the resulting finite sum.

Set Z=epsilon^T P epsilon for six independent signs. The32 representatives
and64 signs give the same average for every function of Z because of
global sign symmetry. BL gives

 int_E exp(-psi)<=16^(-r)exp(E[Z log Z]/2).

Expand Z=trP+2sum_{i<j}P_ij epsilon_i epsilon_j. Distinct off-diagonal
terms have zero product expectation, including those sharing one index;
diagonal terms are constants. Consequently EZ=trP and
EZ^2=(trP)^2+4sum_{i<j}P_ij^2. For an orthogonal projection this is
r(r+2)-2sum P_ii^2<=r(2+2r/3), by sum P_ii^2>=r^2/6. This is exact
finite-sign algebra. No Gaussian replacement or sample extrapolation
occurs. The saved zero-column witness also has an analytic explanation:
three orthogonal pair differences yield Z=2 Binomial(3,1/2), hence
representative multiplicities4,12,12,4 at0,2,4,6.

Under the size-biased probability with mass density Z/r on positive
atoms, Jensen gives E[Z log Z]/r<=log(EZ^2/r)<=log M_r,
M_r=2+2r/3. Its probability mass is1 even when some Z atoms vanish.
Combining constants proves for every original-law rank-r projection

 Q_E(0)<=(2+2r/3)^(r/2)/(32pi)^r, 1<=r<=6.

Absolute integrability is independently available from psi>=sqrt32|xi|
on every subspace. The product characteristic function is that of the
projection of the full original sum of independent Cauchy variables;
Fourier inversion therefore gives the continuous density of that law.
Cos<=1 against exp(-t psi) gives q_E,t(x)<=q_E,t(0), and degree-one
homogeneity gives q_E,t(0)=t^(-r)Q_E(0). Positivity comes from the law,
not from merely taking the real part of a Fourier integral. An arbitrary
deleted-column law does not inherit the frame or moment identities.

## True normal trace and all104 P+B coefficients

For r=3 the bound is Q_new=1/(4096pi^3). Accepted098 proves the pointwise
identity H(z)=int_0^infinity t^2 q_4,t(z)dt=Q_{e-perp}(0)/(2|z|), for
z!=0, using stable homogeneity and the actual continuous line marginal.
This new audit does not infer that trace identity from a peak inequality.
It inserts the new bound into that accepted identity to get
H(z)<=1/(8192pi^3|z|). This holds for every genuine original4D quotient
and every direction of its3D line marginal, not only body-triple strata.
The E098 proof/hash and K104 accepted scientific review are dependencies.

I read104 proof(1)-(12), including restoration before time integration.
Restoration returns the retained scalar Cauchy convolution to the original
full4D quotient. In(6), every complement summand is a positive geometric
factor times the same uniform3D peak. The large P branch is also linear
in that peak: its coefficient is proportional to1/(pi^2 c^3). In(9) and
(12), D and every P coefficient have exactly that common dependence;
all norms' remaining radial, u, Jacobian and cone factors are independent
of the peak. Hence replacing Q_old=1/[pi^2(32/sqrt3)^3] by Q_new multiplies
EVERY P and B norm by lambda=8/(3sqrt3 pi). Their Minkowski sum also
scales by lambda, and each q=6/5 integral allocation scales by
lambda^(6/5). It would be invalid to change B only, or to apply this
constant to the deleted law before restoration. Neither happens here.

The saved104 dyadics, not displayed decimals, are the only budget inputs.
PaidCEH is fixed. The signed remaining allocation is target-paid-newI,
so its lower endpoint subtracts the UPPER endpoints and conversely.
The core allowance is positive conditionally on g0=1/65000, with other
regions unpaid. The complete31-middle constructed upper allocation still
has negative S_plan even if far is set0. This rejects the competitiveness
of that allocation expression; an upper budget exceeding target is not
a lower bound for true I, true model cost or C6. K104's source event
completeness, all930 pairs, gap geometry, face boundaries and Jacobians
are accepted dependencies, not independently regenerated here.

## Rank6 peak and exact scope

For P=I all z_v=6 and the size-biased moment step has no loss. Direct
BL gives integral exp(-psi)<=(sqrt6/16)^6, so the actual full6 peak is
bounded by216/(2^30 pi^6)t^-6. Relative to the saved001 determinant160
upper bound the exact factor is135/4194304. The original full6 density
need not have an isotropic Cauchy spatial profile: the result supplies
no x decay and no far-region budget. The all-r theorem is valid even
where its constants are not sharp. It does not alter the target centered
cube operator into a strong maximal operator or a composition.

The theorem is stronger than the universal sharp-Khintchine coefficient:
their ratio is2sqrt2/pi<1, using pi>3>sqrt8. The author makes no claim
about every possible105 directional cap refinement; none is needed here.
General L1/all positive radii remain unclosed at global g/full I/far/S.
The actual upper C6<=64 is unchanged and no lower-project result is made.

## Source and code qualification

The local special-case proof is self-contained. Public BL references in
the author manifest are bibliographic attribution only; their missing
byte hashes and partial-read status are explicit and not concealed as
fully checked external premises. The heat proof above supplies the exact
case required. The author's run_prepared.py and unexecuted saved checker
were read, not imported/executed. Their three configurations check finite
algebra and Arb arithmetic; those runs do not substitute for this analytic
review. Our one Fraction-only saved checker independently certifies coarse
root/rescaling bounds and signs, preserving all fine saved bytes by hashes.
