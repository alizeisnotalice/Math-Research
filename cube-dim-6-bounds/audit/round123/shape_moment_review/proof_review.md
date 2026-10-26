# E123 independent analytic review

The full-boundary moment mechanism and the analytic Gamma upper baseline
are correct. This review certifies the mathematical interface and the
new baseline scalar arithmetic; it does not audit or repeat the8192-shape
sampling calculation or turn its empirical61.29 into a rigorous bound.
The original body32, unit weights/common Cauchy time and Lebesgue/Fourier
normalization are unchanged.106/027 and accepted original-generator
source results remain explicit pinned dependencies.

## Same-time mixture and complete boundary norm

For positive iid S with Laplace exp(-sqrt(lambda)), Sigma=sum S_i v_i v_i^T
is SPD almost surely. The covariance2t^2Sigma Gaussian mixture has
characteristic function exp(-t psi), so it is the original Cauchy law.
Its positive density integration over t gives C=Gamma(5/2)/(4pi^3)
=3/(16pi^(5/2)) and G=C E[d Q^-5/2]. Tonelli allows infinite values;
the conditional Gaussian generator is never substituted into the source
theorem. Mixture pointwise continuity follows as in106 from finite
negative moments of S and det(Sigma)^(-1/2)<=32^-3 sum S_i^-3.

For EACH fixed SPD Sigma, cube radial coordinates have volume factor
r^5(x.n)drdS=r^5drdS on every unit cube face. Since
int_0^infinity r^5 exp(-r^2Q)dr=Q^-3, diagonal Gaussian integration gives
int_boundary Q^-3 dS=pi^3 sqrt(det Sigma). No single-shape permutation
or sign symmetry is used. The q=6/5 norm of Q^-5/2 is therefore
pi^(5/2)det(Sigma)^(5/12). Multiplying d leaves det(Sigma)^(-1/12).

To prove finite I without assuming it, restrict expectation to the nested
events E_n={n^-1I<=Sigma<=nI}. For fixed n, density/norm integrands are
bounded on the compact boundary and Minkowski applies directly. The
events increase to all SPD matrices. Positive MCT for G_n^q and for
the moment gives ||G||_q(boundary)<=(3/16)M, where
M=E det(Sigma)^(-1/12). Finiteness of M is independently proved below.
This includes all times and all strata; infinite G values on S2 or any
other null set are not discarded to create an estimate.

ONLY the original expectation G is invariant under all coordinate signs
and permutations. Its whole boundary q-integral is384I, giving
I<=((3/16)M)^(6/5)/384. A sampled shape generally lacks that symmetry;
applying384 to one shape's positive-face integral would be invalid.

## Scalar shape-weight optimality

For a positive scalar shape weight w with the needed finite moments,
weighted Jensen followed by the complete boundary identity gives the
factor(Ew)^(q-1)E[w^(1-q)det(Sigma)^(-1/10)]. Hölder bounds this below
by(E det(Sigma)^(-1/12))^q. Equality is attained by w proportional to
det(Sigma)^(-1/12), because its mixed integrand is then the same moment.
Thus the moment budget is optimal in this specific shape-only scalar
weight family. This is not optimality over spatially dependent weights,
all estimates or the actual I/model. For old weight d, strict concave
Jensen gives M^q<(Ed)^(1/5) since d is integrable and nonconstant under
the actual positive iid law; varying one scale on an open event changes
the determinant. It need not imply either budget is competitive.

## Radial factor and self-contained bounded-shape proof

Y_i=1/(4S_i) are iid Gamma(1/2,1). Their sum R is Gamma(16,1), independent
of P~Dirichlet(1/2,...,1/2). Sigma=A(P)/(4R), A=sum vv^T/P_i. Hence
det(Sigma)^(-1/12)=2sqrt(R)Z and
M=D E Z, D=2Gamma(33/2)/Gamma16. The half-integer factorial identity
reduces D to(300540195/67108864)sqrt(pi). No random radial sample is needed.

Cauchy-Binet gives det(sum t_i v_i v_i^T)=sum_J c_J prod_(i in J)t_i,
c_J=det(V_J)^2>=0. At t_i=1, sum c_J=32^6. Differentiating there via
Jacobi's determinant formula gives normalized inclusion probability
v_i^T(32I)^-1v_i=6/32=3/16. This is exact leverage, not an enumeration
of basis subsets. Omit zero-weight subsets when applying weighted AMGM;
the remaining weights sum1. Consequently

 det A>=32^6 prod_i P_i^(-3/16)>=32^12,

where the last step uses prod P_i<=32^-32. At uniform P, A=1024I,
so Z=det(A)^(-1/12)<=1/32 and equality is attained there. Interior
Dirichlet shapes have positive coordinates almost surely. Thus Z is
positive/bounded, M<=D/32 is finite and the Minkowski argument is closed
without a circular assumption on I or a numerical p1 computation.

Retaining the product before that uniform bound gives
Z<=32^-1/2 prod P_i^(1/64). The exact Dirichlet beta integral has total
exponent32/64=1/2, so

 E prod P_i^(1/64)=Gamma16/Gamma(33/2)
                  *[Gamma(33/64)/Gamma(1/2)]^32.

Multiplying D cancels its radial Gamma ratio and proves the genuine
complete analytic expectation upper M<=M_Gamma,
M_Gamma=(2/sqrt32)[Gamma(33/64)/Gamma(1/2)]^32.
This bound is distinct from any MC mean/confidence assertion. An interval
evaluation of this explicit constant certifies the analytic baseline,
not the actual moment's location between its two numerical endpoints.

## Threshold and full same-generator positive-part mass

The complete moment coefficient obeys I_M^5=(3/2^59)M^6. For globally
certified g_L, the old same-G positive-part transfer gives
U_M=6^(6/5)I_M/g_L^(6/5)=(9/8)^(6/5)(M/g_L)^(6/5)/384.
At the UNPROVED g0=1/65000, U_M<64 is equivalent to

 M^6<(8/9)^6*24576^5/65000^6
    =2^83/(3^7*65000^6)
    =36893488147419103232/629200415790081024169921875.

Our exact check covers all equivalent constants; main MATLAB's three
previous symbolic identities are provenance, not expectation evidence.
The new scalar configuration evaluates M_Gamma and this threshold with
outward endpoints to decide only the analytic baseline's competitiveness.

For c>0, G has degree-5 and nonnegative cube-radial integration gives
int(G-c)_+=(5/6)c^-1/5*384I=320c^-1/5I. The new complete I estimate
therefore proves global L1 mass of the actual spatial positive part,
including arbitrary short/long times and all angular singularities.
It does not truncate conditional Gaussian shapes before mixing. With a
certified0<c<g_L, normalization by64(g_L-c) majorizes the cube and has
mass5I c^-1/5/(g_L-c). Its log derivative is-1/(5c)+1/(g_L-c), uniquely
zero at c=g_L/6, with the endpoints divergent and the stated minimum.
The established whole-space positive-source/obstacle argument for this
SAME Cauchy generator then gives all generalL1 inputs and all positive
common cube radii; no conditional Gaussian or atomic-input restriction
is introduced. A competitive matrix moment upper and globalg0 are still
separate missing certificates. A coarse Gamma I upper is complete but
may fail competitiveness without excluding actual M, I or the model.

## Evidence qualification

The author MC output is explicitly diagnostic and not independently
averaged, resampled or statistically certified here. Its source/input/
proof/release/started/receipt/chunks are byte-preserved only. The claimed
8192 samples, empirical61.29 and lack of fallbacks remain historical
author evidence, not part of this analytic acceptance. Our one new
scalar run is at most one CPU and calls no author/old kernel program,
MC, native, optimization, MATLAB, Astra, recursive agent or lower project.
Parent retains final proof responsibility.
