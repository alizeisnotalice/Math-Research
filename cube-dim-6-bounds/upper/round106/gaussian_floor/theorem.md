# R106 Gaussian floor: exact representation and an order-of-minimization obstruction

Status: author derivation read and released by the primary; the one supporting verification ran once and passed. Independent audit remains for the primary to arrange. This concerns the existing body32 upper-bound model only. No lower-bound research and no new C6 bound are claimed.

## 1. Scope, normalization, and prior-use check

Let V contain all 32 vectors (1,s2,...,s6), sj in {−1,1}, each with weight one. The original Cauchy semigroup is defined by characteristic function exp(−t psi(xi)), psi(xi)=sum_v |v dot xi|, and G=integral_0^infinity p_t dt. Its homogeneity is degree −5, with the existing q=6/5 and upper-bound accounting unchanged.

The independent clone was attempted at the new, previously absent path Math-Research-cube6-sol-gaussian106. The local clone operation exited 128 after about 1.46 seconds: its independent .git was created, but checkout had missing promisor blobs. The available project checkout is used transparently; it is not called a complete checkout. No stage/commit/push was performed by this agent. Initial clone status already displayed index deletions caused by failed checkout. All new files belong only to upper/round106/gaussian_floor.

Limited duplicate check: upper/round006/alternative/report.md section 6 already gives the general alpha positive Gaussian-subordinator representation and an alpha=3/2 sampling pilot. The same file explicitly stops that pilot and does not certify a cube floor. Searches of upper/theory, rounds001–006,091,093,100,101 for Gaussian/subordinator/determinant/Dirichlet found no use of the all-cube random-matrix minimum below. Thus the positive representation itself is not claimed new; the obstruction and sign-orbit alternative are the present contribution.

## 2. Exact positive representation (including the factor 1/4)

Take independent positive variables S_v with Laplace transform E exp(−lambda S_v)=exp(−sqrt(lambda)). One explicit law is

f_S(s)=(2 sqrt(pi))^(−1) s^(−3/2) exp(−1/(4s)), s>0;
equivalently S=1/(4Y), Y~Gamma(1/2,1), or S=1/(2Z^2), Z standard normal.

Set Sigma=sum_v S_v vv^T. Since sum_v vv^T=32 I and all scales are positive, Sigma is positive definite almost surely. For each fixed t>0, conditionally Gaussian with covariance 2t^2 Sigma has characteristic function exp(−t^2 xi^T Sigma xi); averaging gives exactly exp(−t psi(xi)). This uses the same t for every direction, and **t^2 Sigma**, not t Sigma. A fixed sampled Sigma by itself is not the original semigroup.

Its conditional density is

p_t(x|Sigma)=(4 pi t^2)^(−3) det(Sigma)^(−1/2) exp(−Q_Sigma(x)/(4t^2)),
Q_Sigma(x)=x^T Sigma^(−1)x.

The mixture density equals the original p_t at every x. For example, with m=min_v S_v, Sigma>=32m I gives det(Sigma)^(−1/2)<=(32m)^(−3)<=32^(−3) sum_v S_v^(−3), which has finite expectation. This supplies domination and continuity for the mixture density. Equality follows from uniqueness of the characteristic function and continuity.

Nonnegative Tonelli gives, allowing +infinity,

G(x)=C E[det(Sigma)^(−1/2) Q_Sigma(x)^(−5/2)],
C=Gamma(5/2)/(4 pi^3)=3/(16 pi^(5/2)).                 (1)

Indeed u=Q/(4t^2) yields t^(−6)dt=−16 Q^(−5/2)u^(3/2)du, so the time integral is 16 Gamma(5/2)Q^(−5/2); multiplication by (4pi)^(−3) gives (1). The proposed denominator 2pi^3 is a factor two too small. The denominator 4pi^3 also agrees with round006's general-alpha formula at alpha=1.

At x=0 both sides are infinite. On the actual S2 union of 496 pairplanes, the existing infinite-G result is compatible with (1): positive expectation may diverge even though every sampled SPD integrand is finite at nonzero x. No finite sample value or matrix event makes G finite on those strata.

The original Cauchy generator and source conditions remain those of the original exact G. Individual conditional Gaussian generators cannot be substituted into those conditions, and truncation must not be interchanged with mixing. This note supplies only lower estimates for the same G.

## 3. The legal whole-cube floor and its strict ceiling

For SPD Sigma let

M(Sigma)=max_{c in {−1,1}^6} c^T Sigma^(−1)c,
H(Sigma)=det(Sigma)^(−1/2) M(Sigma)^(−5/2).

A convex quadratic reaches its maximum on the cube at a corner; a direct proof writes x as a convex combination of corners and applies convexity. Thus for every x in [−1,1]^6, the integrand in (1) is at least H(Sigma). Therefore

B_Gauss:=C E H(Sigma) <= inf_cube G <= essinf_cube G=g. (2)

The expectation is strictly positive because H>0 almost surely. It is finite by the estimate below. This proves legality on the entire cube, including singular strata, and hence for the essential floor. Matrix-event lower bounds that discard the complement are further lower bounds on B_Gauss.

**Proposition (intrinsic ceiling for this order of operations).** Put T=tr(Sigma)/6=sum_v S_v. For every SPD Sigma,

H(Sigma) <= 6^(−5/2) T^(−1/2).                        (3)

Proof. Averaging corner quadratics gives tr(Sigma^(−1)), so M>=tr(Sigma^(−1)). Let a1,...,a6 be the eigenvalues of Sigma^(−1), e1=sum a_i, e5=sum_i product_{j!=i}a_j. Then

[det(Sigma)^(−1/2) (tr(Sigma^(−1)))^(−5/2)]^2
 = product_i a_i / e1^5.

Also tr(Sigma)=e5/product_i a_i. The elementary Maclaurin inequality e5<=6(e1/6)^5=e1^5/6^4 implies product a/e1^5<=1/(6^4 tr(Sigma))=1/(6^5 T), proving (3). All steps preserve inequality direction.

For completeness, the particular Maclaurin inequality used here follows by pair averaging at fixed sum: with four remaining nonnegative entries r, e5(a,b,r)=(a+b)e4(r)+ab e3(r). Replacing a,b by their common average increases this expression. A maximum on the fixed-sum compact simplex cannot have fewer than five positive entries (then e5=0), and cannot have exactly five (averaging its zero with a positive entry strictly increases e5). An interior maximum must have every pair equal, so its value is 6(e1/6)^5. This proves precisely the needed inequality without any matrix probability assumption.

The 32 independent S_v have T distributed as 32^2 S: its Laplace transform is exp(−32 sqrt(lambda)). The negative moment E S^(−1/2)=2/sqrt(pi), obtained from the Laplace transform by the positive Gamma-integral, gives E T^(−1/2)=1/(16sqrt(pi)). Consequently

0 < B_Gauss <= 3/(256 pi^3 6^(5/2))
               < 1/165888 < 1/65000.                  (4)

The first strict rational upper estimate uses only pi>3 and sqrt(6)>2. In particular, no exact expectation integration, shape conditioning, symmetry of the distribution, or positive matrix event can enable the particular lower bound B_Gauss to certify g0=1/65000. This is stronger than a failure of an all-scales-in-one-box event. It does not rule out the original G, a positive S with other inputs, or a different order of minimization and averaging. It is not a lower bound for C6.

A tempting determinant-weighted Jensen collapse also cannot repair this loss. If m0=E det(Sigma)^(−1/2) and m1=E[det(Sigma)^(−1/2)tr(Sigma^(−1))], symmetry yields the valid floor C m0^(7/2)/m1^(5/2). Jensen bounds this by C E[det(Sigma)^(−1/2)(tr(Sigma^(−1)))^(−5/2)], which obeys exactly the same ceiling (4). Thus these two scalar moments are not the proposed remedy.

## 4. Exact radial self-normalization

Write Y_v=1/(4S_v), R=sum_v Y_v, and P_v=Y_v/R. Then R~Gamma(16,1), P~Dirichlet(1/2,...,1/2) on the 31-simplex, independently. With A(P)=sum_v P_v^(−1)vv^T one has Sigma=A(P)/(4R). Homogeneity gives

G(x)=2C Gamma(33/2)/Gamma(16)
     * E_P[det(A(P))^(−1/2)(x^T A(P)^(−1)x)^(−5/2)]. (5)

The same formula holds with max-corner quadratic for B_Gauss. This analytically removes one radial variable and all its heavy tail; it does not turn the remaining shape integral into low dimension. Boundary singularities of the Dirichlet integral still encode the actual singular strata. In particular (5) cannot evade (4) if its integrand has already been minimized over the whole cube.

## 5. A concrete different matrix control: reflect before minimizing

For each diagonal sign matrix D_epsilon, epsilon in {−1,1}^6, define the finite orbit average

F_Sigma(x)=2^(−6) sum_epsilon
   det(Sigma)^(−1/2) [x^T D_epsilon Sigma^(−1)D_epsilon x]^(−5/2),
H_orbit(Sigma)=inf_{x in cube, x!=0} F_Sigma(x).

The law of Sigma is invariant under conjugation by all D_epsilon (body directions are merely permuted up to sign). Hence E F_Sigma(x) is precisely the integrand expectation in (1). Therefore

B_orbit=C E H_orbit(Sigma) <= g.                       (6)

The infimum is an actual minimum on the cube boundary: F is homogeneous of degree −5, continuous off 0, and tends uniformly to infinity as x tends to 0. It is a measurable continuous function of SPD Sigma. Its expectation is finite by comparison with G at any one regular point outside S2, and is positive. The minimizer of the sign average need not be a corner; one must not reapply the single-quadratic corner argument to F.

This quantity preserves a common x while the 64 reflected covariances are averaged. It satisfies H_orbit>=H but, unlike H and its trace/Jensen relaxation, it need not obey (3). Here is one finite, directly verifiable example, with no random sampling.

Take v=(1,1,1,1,1,1), Sigma=I+100vv^T. This belongs to the actual positive body32 cone: all scales equal 1/32, with an additional 100 on the v direction. Then det(Sigma)=601, T=101, and Sigma^(−1)=I−(100/601)vv^T. For rho=100/601,

F_Sigma(x)=601^(−1/2) * 2^(−6) sum_epsilon
     (|x|^2−rho(epsilon dot x)^2)^(−5/2).              (7)

The minimum in (7) reduces to five one-dimensional closed branches. To justify this, fix r^2=sum_i x_i^2 and write u_i=x_i^2. In the uniformly convergent positive expansion of (1−rho Z^2/r^2)^(−5/2), Z=sum epsilon_i x_i, every even Rademacher moment decreases when two u_i with fixed sum are made more unequal. Indeed for their pair sum W,

E W^(2l)=sum_{j=0}^{floor(l/2)} binom(l,2j)
             (u_i+u_j)^(l−2j) 4^j (u_i u_j)^j.

Each coefficient is nonnegative; averaging the independent remaining signs keeps all coefficients nonnegative. Successive transfers therefore minimize all the moments at a packed vector: k coordinates equal one, one coordinate y in [0,1], the others zero. Homogeneity confines the minimum to max |x_i|=1, so k=1,...,5 suffices (k=5,y=1 includes the all-ones vector).

The resulting pre-determinant branch is

A_k(y)=2^(−k−1) sum_{ell=0}^k binom(k,ell)
 [ (k+y^2−rho(k−2ell+y)^2)^(−5/2)
  +(k+y^2−rho(k−2ell−y)^2)^(−5/2) ].                  (8)

The verification source tests the concrete sufficient statement

A_k(y)>1/20 for all five branches and all y in [0,1]. (9)

It partitions each y interval into 64 closed rational subintervals. Every denominator is a convex quadratic in y, so its exact rational upper bound on a subinterval is the larger endpoint value. Substituting those upper bounds into the decreasing power gives a lower bound for each summand, and 192-bit Arb then certifies the finite sums. No interval enclosure assumes a derivative sign.

The sole released execution passed all 320 closed subintervals at 192-bit Arb precision, proving (9) and therefore

H_orbit(Sigma)>1/(20sqrt(601))
                >6^(−5/2)/sqrt(101)>=H(Sigma),       (10)

where the middle comparison is the exact integer inequality 400*601 < 6^5*101. Thus the alternative is a concrete control with a strict deterministic improvement on an actual body32 covariance. Continuity also gives an open covariance neighborhood with that strict gap. This does not certify E H_orbit or g0; it identifies a mechanism that escapes the proved ceiling rather than resampling the same blocked whole-cube minimum.

For this example the original H is exactly 1/(sqrt(601)6^(5/2)): a balanced corner with three plus and three minus signs has v dot c=0, so the maximum corner quadratic equals 6. Thus (9) improves that single-matrix floor by a factor greater than four (6^5>80^2), while also exceeding its looser trace ceiling in (3).

## 6. What enters the research ledger

The actual Cauchy Gaussian representation has the coefficient in (1). The whole-cube samplewise minimum is a legal positive floor but has the strict structural ceiling (4), which excludes g0 for that particular mechanism. The radial elimination (5) is exact. The orbit-before-minimum control (6) is a different legal floor, with a finite deterministic test (8)–(10). Random expectation certification for (6), a competitive complete I, and positive S are all unpaid. Formal C6 upper bound stays 64, and the stopped lower-bound research remains untouched.

## 7. Actual verification and reproduction

Primary release: primary_release.json, SHA256 77b9c0cc9ed9ad68ed32edcf67a0b80d9f4885387d711da29ff88ca1a9e8d4e0. It binds input SHA256 18c9ea9b6dc36491a3ae25625213f4b18fbda847d1da08b96408b4a01a44bd37 and source SHA256 27ccca5717a0feaf12e3f4dd4dc60b5f678fee317c467a0cee0d6a1582db93e5. The primary's static review covered the inverse-eigenvalue inequality, positive-moment packing, uniform convergence, actual matrix-cone membership, and finite endpoint bounds before execution.

The one actual execution exited 0 (tool chunk dcdbb1), with program receipt wall time 0.374116291 seconds, zero child processes, no persistent session, and all 320 checks passed. The source's SIG_DFL 177-second alarm remained armed through receipt and exit; it was retained for this finite short-feedback batch, not described as a current user runtime ceiling. The output contains full finite interval rows and hashes. The smallest displayed row lower bound was about 0.09548241606, far above the accepted exact 1/20. Selecting that display row used a float and had no role in acceptance; every row was accepted by an Arb comparison.

The ceiling in (4) was enclosed by Arb at approximately 4.286013605856154e−6. The sufficient orbit example bound 1/(20sqrt(601)) is approximately 0.00203954254112, while its original trace ceiling is approximately 0.00112839508848. These are supporting interval values; the exact statements are (4), (9), and the positive integer comparison in (10).

Reproduction, from the independent clone root:

```
/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python cube-dim-6-bounds/upper/round106/gaussian_floor/check_gaussian_floor.py
```

This command has already been run once under the recorded release. Further executions require the primary's scope decision; the present deliverable does not authorize another run by implication. No Monte Carlo, old kernel point, 091 baseline, MATLAB, lower-bound work, staging, commit, or push was performed.
