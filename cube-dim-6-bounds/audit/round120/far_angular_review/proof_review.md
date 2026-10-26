# E120 independent analytic review

PASS for the stated positive-mixture interface, specified Jensen-budget
obstruction and complete coarse actual-far upper bound. None supplies a
lower bound on actual I_F or excludes the original potential/model/C6.
I read the full120 proof, producer, observer and saved checker;106's
positive representation/normalization,027's same-generator mechanism
and boundary factor,059's actual region, plus the complete saved111
comparison. Prior108 peak and exact stored pi/111 account are dependencies,
not rerun calculations. Historical037/098/104/107/006/003 references mark
prior-use boundaries; none supplies a new unproved step of the120 proof.

## Positive original-law mixture and time integral

Take independent S_v>0 with Laplace transform exp(-sqrt(lambda)), and
Sigma=sum S_v vv^T. The original32 frame has sum vv^T=32I, so Sigma is
SPD almost surely. E exp(-t^2 xi^T Sigma xi)=exp(-t sum|v.xi|) proves
the conditional covariance is2t^2 Sigma, and its mixture is the SAME
common-time Cauchy law. It is not a conditional Gaussian generator
substituted into the original source/obstacle argument.

For m=min S_v, d=det(Sigma)^(-1/2)<=(32m)^(-3)
<=32^(-3)sum S_v^(-3). The explicit S=1/(4Y), Y~Gamma(1/2,1), gives
finite negative third moment. Thus m0=E d is finite and positive; the
mixture densities have an integrable uniform dominating factor and are
continuous. Characteristic-function equality upgrades their a.e. equality
to equality everywhere. At time1 and x=0,
a=p1(0)=(4pi)^(-3)m0, so m0=(4pi)^3a exactly.

Nonnegative Tonelli integrates all t first. For Q=x^T Sigma^(-1)x,
the substitution u=Q/(4t^2) gives int t^-6 exp(-Q/(4t^2))dt
=16 Gamma(5/2)Q^(-5/2); multiplying(4pi)^(-3) gives
C=3/(16pi^(5/2)). At nonzero x each SPD integrand is finite; the
expectation may be infinite and is not discarded. The identity is
extended nonnegative and introduces no oscillatory interchange.

Under probability d dP/m0, Jensen with q=6/5 gives
G(x)^q<=C^q m0^(q-1) E[d Q^(-5q/2)]
=C^q m0^(1/5) E[d Q^-3]. Apply Jensen to bounded truncations of
Q^-5/2 and then use monotone convergence if needed. This is after the
exact time integral; no per-time Lq or time-Minkowski step is hidden.

## Critical-degree spatial identity, including normal measures

Let A_F={z!=0:dist(z,S2)>=rho||z||infinity}. Since S2 is a finite union
of LINEAR actual body-pair planes, A_F is a measurable cone, exactly the
radial extension of all actual far faces. On the positive face it agrees
with059's F by the31 unsigned nearest-distance identity; in signed space
the definition continues to use all496 actual planes. Distances, plane
intersections, tube overlaps and all >=rho boundary conventions are kept.

For a cube face, write z=r x, where x lies on its open flat face. The
volume Jacobian is r^5(x.n)dS(x)dr=r^5dS(x)dr because x.n=1 on every
face of[-1,1]^6. No extra factor6 appears. Equivalently the radial map
theta=x/|x| satisfies dS_sphere=(x.n)|x|^-6 dS_cube=|x|^-6dS_cube.
For f=Q^-3, f(theta)=|x|^6f(x), so the two factors cancel EXACTLY.
Consequently integral_(cube boundary intersect A_F)f dS equals
integral_(S5 intersect A_F)f dS. The same identity follows by comparing
the two radial decompositions of a critical homogeneous annular integral;
each radial integral contributes log L. Face edges/ties have surface
measure0 and need no additional normal measure.

A centered Gaussian with covariance2Sigma has density
(4pi)^(-3)d exp(-Q/4). Its radial integral is
int_0^infinity r^5 exp(-Q(theta)r^2/4)dr=64/Q(theta)^3.
Therefore its angular density is d Q(theta)^-3/pi^3 on S5. This integrates
to1 and uses surface area|S5|=pi^3, not a probability-normalized sphere
measure inserted a second time. Thus for every fixed SPD Sigma,

 int_(S5 intersect A_F)Q^-3 dS=pi^3 sqrt(det Sigma) nu_Sigma(A_F).

The law of Sigma is invariant under coordinate reflections/permutations:
they permute the32 rank-one matrices vv^T, and the S_v are iid. A_F has
the same invariance. Hence the EXPECTED integrand has384 congruent far
face pieces:6 choices of fixed coordinate,2 signs,32 free orthants.
The identity does not assert that one sampled ellipsoid itself has these
symmetries. Tonelli first justifies averaging the full boundary identity;
the factor d then cancels sqrt(det Sigma). The conditional time1 Gaussian
mixture is the original time1 law, so E nu_Sigma(A_F)=Pr(X1 in A_F)=P_F.
The time1 law is absolutely continuous; X1=0 and maximal-coordinate ties
have probability0. No singular angular set is removed to manufacture a
finite integral: the full integral is controlled by this probability.

This proves actual I_F<=J_F with
J_F=C^q m0^(1/5)pi^3 P_F/384. Raising to5 and substituting m0 gives

 J_F^5=(3pi^3 a/2^53)P_F^5.

The fifth power constant follows from C^6*(pi^3/384)^5*(4pi)^3a;
pi exponents are-15+15+3=3, powers of2 are-24-35+6=-53,
and powers of3 reduce from6 to1. This is an integral budget for the
complete actual far region and all time, not a pointwise isotropic
Cauchy heat shape or an independently changed source theorem.

## Actual far probability at the single fixed R

Each coordinate marginal of the original X1 is Cauchy(32), irrespective
of coordinate dependence. The six-coordinate union bound is
Pr(||X1||infinity>R)<=12 atan(32/R)/pi<384/(pi R).
At R=256 this is3/(2pi). If X1 is outside A_F and has maxnorm<=R,
some actual unordered pair plane has orthogonal projection of norm<rho R.
There are496=32*31/2 such pairs; every pair has rank2 and its actual
original projected law has dimension4.108 bounds its peak by
(14/3)^2/(32pi)^4. Its ball event is bounded by this peak times
(pi^2/2)(rho R)^4. Summing all496 events only uses a union upper bound;
no overlap or independence assumption is made. It is crucial not to use
31 unsigned planes globally or992 ordered labels as distinct events.
At rho R=4 the exact coefficient is1519/(1152pi^2).

Thus P_F>1-3/(2pi)-1519/(1152pi^2)>3665/10368=P0, since the expression
is increasing in pi>0 and pi>3. The same fixed-R refined lower bound
uses the complete pi LOWER endpoint; it is not another radius experiment.

## Envelope obstruction, actual coarse upper and account directions

Cauchy-Schwarz and the exact frame imply psi(xi)<=sqrt32 sqrt(32|xi|^2)
=32|xi|. Positive Fourier inversion AT0 gives
a>=(2pi)^(-6)|S5| Gamma(6)/32^6=15/(2^33pi^3).
This is an analytic coefficient witness inside upper-model budget
analysis, not a numerical peak refresh or a C6/lower-project theorem.
Together with the strict probability lower bound,

 J_F^5>(45/2^86)(3665/10368)^5=W5.

The saved exact W5-R_U^5 is positive, using the full111 ref108 endpoints;
both endpoint comparisons are checked. Therefore this SPECIFIC
det-weighted Jensen upper-budget expression cannot fit the current
allocation even with exact a and exact P_F. It does not imply I_F>=W,
does not refute actual G/model/C6, and does not exclude a bound preserving
shape compensation before Jensen. The loss is exactly localized at that
Jensen step. The residual target remains conditional on globally
unproved g0=1/65000.

Conversely108 gives a<=216/(2^30pi^6) and P_F<=1, so
actual I_F<=J_F<=(81/(2^80pi^3))^(1/5). This is a valid complete,
noncompetitive actual-far upper bound. Its two saved endpoints enclose
the UPPER-BOUND COEFFICIENT, not actual I_F from below. Comparing this
coefficient with111's allowance cannot establish actual-I exclusion.
Our saved check verifies allowance copies from final/ref108, target
subtraction directions and exact powers without rerunning111 geometry.

## Execution and final scope

Author evidence records one science run with two completed configurations
and one saved-only check, not a continuation. The observer covers the
science child's spawn through final receipt, exit and reap; its own final
receipt/exit are outside that scope. All17 source pins and final bytes
remain bound. Our independent checker does not import these generators,
evaluate Gaussian samples, old functions, MC, native, kernels or grids.

027's whole-space spatial positive-part source and generalL1/all common
positive radii remain the SAME original generator/G. Gaussian conditioning
is used only to estimate G and its face integral; truncation/source
operations are not interchanged with mixing. Competitive globalg0, far
allocation, complete positive S and a new Cn remain unpaid. No lower
project work is undertaken. Parent retains final proof responsibility.
