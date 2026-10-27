# Seven-dimensional finite-source interface and singular strata

2026-10-27. Round001 theory. This is an independent dimension-seven derivation, using the read-only second-review files listed in `sources.json` as proof provenance. The source statuses are not premises. Every finite computational statement below was actually checked by `check.py`; every assertion about the continuum has its analytical argument here. No MATLAB was invoked by this agent, no old three-/four-dimensional Green nodes were recomputed, and no g/I search was performed.

**Status.** The true semigroup, obstacle, whole-space finite source, one-domain/all-positive-scales/general-L1 comparison, exact singular set, continuous finite Green function off that set, whole-face short-time L^(7/6) budget, and conditional g/I cost are established. There is no new numerical C7 upper bound: neither a useful certified global g0 nor a certified full I0 has been supplied. The finite structural computation is auxiliary evidence for counts and constants, not a substitute for the analytical proof.

## 1. Fixed model, genuine density and domains

All measures are Lebesgue unless indicated. Fix

    V={(1,s2,...,s7): sj in {-1,1}}, m=64,
    psi(xi)=sum_{v in V}|v·xi|, L=-psi(D),
    X_t=sum_{v in V} v Z_v(t),

where the 64 variables are independent centred one-dimensional Cauchy variables with the **same** scale t and density f_t(z)=t/[pi(t²+z²)]. Vectors are unnormalised and every weight is exactly1. Its characteristic function is exp(-t psi), so X_t yields a positive, conservative, symmetric convolution semigroup P_t. This is one process, not a product of maximal operators with separate times.

The exact Gram identity is VV^T=64 I7. Consequently

    8|xi| <= psi(xi) <= 64 sqrt(7)|xi|.

The lower estimate follows from l1>=l2 on the64 scalar products; the upper estimate uses |v|=sqrt(7). Therefore Dom_2(L)=H1 and the form domain is H^(1/2), with inhomogeneous equivalent norms. No H1-to-continuous embedding is asserted.

An explicit basis B with |det B|=512 is obtained from the Sylvester Hadamard matrix H8: remove its first row/column, multiply columns by their first entries, and then multiply rows so its first column becomes all1. The actual integer B appears in `results.json`. All its columns lie in V. For any target v, diag(v)B is another direction basis, contains v as its first column, and has determinant absolute512. The finite check verifies these64 bases exactly. Thus

    q_B,t(x)=|det B|^(-1) product_{j=1}^7 f_t((B^(-1)x)_j)

is strictly positive, continuous, C0 and of mass1. Convolution with the independent remaining57 variables gives the canonical density p_t, with

    p_t>0, p_t in C0, integral p_t=1,
    ||p_t||infinity <= C t^(-7), C=1/(512 pi^7).

Strict positivity follows from q_B,t>0 and a probability convolution. A C0 function convolved with a finite measure stays C0. Its Fourier transform is integrable because psi>=8|xi|, so the density also agrees everywhere with the Fourier inversion convention (2pi)^(-7) integral exp(ix·xi-t psi(xi)) dxi. In particular p_t is smooth in x and jointly continuous for t>0. Scaling of all64 Cauchy variables gives

    p_t(x)=t^(-7)p_1(x/t).

Dilation continuity in L1 gives positive-time L1 continuity of p_t. There is no hidden independent scale or change of direction weight.

For phi smooth compactly supported, the real generator is

    Lphi(x)=(1/(2pi)) sum_v integral_R
      [phi(x+sv)+phi(x-sv)-2phi(x)]/s² ds.

The second difference is bounded by min(s²|v|²||D²phi||infinity,4||phi||infinity). Its L1 bound is obtained with the analogous L1 derivative/translation estimates. Splitting the s integral at a fixed positive number proves Lphi in L1 intersect C0 intersect Linfinity. For a fixed nonnegative smooth radial cutoff chi, equal1 near0 and nonincreasing in radius, write chi_R(x)=chi(x/R). Homogeneity gives

    ||Lchi_R||infinity=R^(-1)||Lchi||infinity,
    ||Lchi_epsilon||1=epsilon^6||Lchi||1.

These estimates and standard Plancherel, weak L2 compactness, local Rellich compactness, Tonelli, monotone/dominated convergence and positive-distribution-to-Radon-measure are the functional-analytic tools used below.

## 2. Obstacle decomposition in dimension7

Let f>=0 be bounded, compactly supported and nonzero, F=||f||infinity, M=||f||1, kappa>0. For 0<t<=1 set

    V_t=t sum_{j>=0}P_{jt}f,
    u_(t,0)=0,
    u_(t,k+1)=(P_t u_(t,k)+t(f-kappa))_+.

The heat bound gives ||P_s f||infinity<=min(F,CMs^(-7)). Hence the j>=1 series converges uniformly for each t, V_t-tf is C0, and positive Tonelli gives P_tV_t=V_t-tf. The iteration is increasing, bounded by V_t, and its limit u_t solves the fixed-point equation. The sharper bound

    0<=u_t<=(V_t-kappa t)_+

follows at every iteration from P_tV_t=V_t-tf. Outside the compact essential support of f, V_t is C0; this bound proves u_t has essential compact support for each fixed t, and is L1 intersect L2. Neither V_t in L1 nor a uniform support radius is needed.

Let L_t=(P_t-I)/t and sigma_t=f+L_tu_t. The positive and zero branches of the fixed point give

    0<=sigma_t<=kappa,
    sigma_t=kappa a.e. on {u_t>0}.

Only after u_t is known to be L1 do we use conservativity to conclude integral sigma_t=M. Pairing f=sigma_t-L_tu_t with u_t gives

    <u_t,-L_tu_t>+kappa||u_t||1=integral f u_t.

The right Riemann sum of the decreasing envelope min(F,CMs^(-7)) yields, uniformly for t<=1,

    ||V_t||infinity <= F+(7/6)(CM)^(1/7)F^(6/7)=:A_f.

It follows that

    ||u_t||1<=MA_f/kappa,
    ||u_t||2²<=MA_f²/kappa,
    ||sigma_t||2²<=kappa M,
    ||L_tu_t||2<=||f||2+sqrt(kappa M)=:A_2.

The L1 control comes from the energy identity, not from an assumed integrable potential.

The multiplier psi_t=(1-exp(-t psi))/t satisfies psi_t>= (1-exp(-1))Q whenever psi>=Q and t<=1/Q. Thus

    integral_{psi>=Q}|u_t hat|²
      <= A_2²/[(1-exp(-1))²Q²].

For each fixed Q the low-frequency projection has a uniform H1 bound. Local Rellich for that projection, followed by the uniformly small high-frequency tail and a diagonal argument, gives a sequence t_j down0 such that u_tj converges strongly in L2_loc and a.e. to u, and weakly in L2; sigma_tj converges weakly in L2 to sigma. Fatou gives u>=0 in L1 intersect Linfinity. For smooth tests, L_tphi->Lphi in L2, so Lu=sigma-f distributionally. Fourier transformation then gives psi uhat in L2, hence u in H1.

Weak closedness of the order interval gives 0<=sigma<=kappa. On D={u>0}, a.e. convergence implies u_tj>0 eventually. For any bounded, bounded-support measurable test supported on D, dominated convergence of sigma_tj to kappa on D and weak L2 convergence imply sigma=kappa there. This is the required saturation argument; weak convergence alone would not imply it.

For chi_R increasing to1, weak convergence first gives integral chi_R sigma<=M and hence sigma in L1 with total mass at mostM. Next the equation and the large cutoff estimate imply

    |<Lu,chi_R>|=|<u,Lchi_R>|
       <=R^(-1)||u||1||Lchi||infinity ->0.

Therefore integral sigma=integral f=M. We have proved

    f=sigma-Lu,
    u>=0 in L1 intersect Linfinity intersect H1,
    0<=sigma<=kappa, sigma=kappa on D,
    integral sigma=M, kappa|D|<=M.

This recovers mass rather than assuming weak limits conserve it. No statement Lu=0 on D^c is made: nonlocal jump-in can be positive there.

## 3. Green function and whole-space finite source

Define the extended nonnegative function G(x)=integral_0^infinity p_t(x)dt. It is even and homogeneous of degree-6. Put U_T=integral_0^T p_tdt and W_T=integral_T^infinity p_tdt. Then

    ||U_T||1=T,
    ||W_T||infinity<=C/(6T^6).

Choosing T=(C/(3lambda))^(1/6) makes the second bound lambda/2, and Markov's inequality gives

    |{G>lambda}|<=2(C/3)^(1/6)lambda^(-7/6).

Thus G is weak-L^(7/6). For each finite-volume set E, direct Tonelli also gives

    integral_E G<=integral_0^infinity min(1,C|E|t^(-7))dt<infinity.

In particular G is locallyL1 and finite a.e.; it need not be pointwise finite on all lower-dimensional sets. For every c>0,

    H_c=(G-c)_+ in L1,

because its integral equals integral_c^infinity |{G>s}|ds. All space and the full unbounded tails have already been included.

For phi smooth compactly supported,

    integral G|Lphi|
      <=integral_0^infinity min(||Lphi||infinity,Ct^(-7)||Lphi||1)dt<infinity.

The semigroup derivative and P_Tphi(0)->0 therefore give LG=-delta0 as a whole-space distribution. Specifically integral G Lphi = integral_0^infinity d(P_tphi(0))/dt dt=-phi(0). Fourier inversion justifies the derivative; psi|phihat| is integrable.

Positive Tonelli and the convolution law give P_tG=integral_t^infinity p_sds, bounded for t>0. Since G=P_tG+U_t a.e., Jensen and the1-Lipschitz positive-part function give

    P_tH_c >= (P_tG-c)_+ >= H_c-U_t.

Test against nonnegative smooth compactly supported phi and divide by t. The left side converges to <H_c,Lphi>, because H_c is L1 and (P_tphi-phi)/t->Lphi uniformly. The latter uniform convergence follows from dominated Fourier inversion with psi|phihat|. The right side converges to -phi(0), since U_t/t converges weakly to delta0. Consequently

    mu_c := LH_c+delta0 >=0

on the **whole space**. This positive distribution is a Radon measure. No direction subspace, free boundary, tail or interface has been removed, and no regularity of the level set {G=c} is required.

The large cutoff gives

    integral chi_R dmu_c=1+integral H_c Lchi_R ->1.

Positive monotone convergence then proves mu_c(R7)=1. Set h_c=min(G,c)=G-H_c as a locally integrable equivalence class. It is bounded, and mu_c=-Lh_c. All pairings are legitimate because Lphi is L1, h_c is bounded, and the G/H_c pairings above converge absolutely. The small cutoff gives

    0<=integral chi_epsilon dmu_c
      <=c epsilon^6||Lchi||1 ->0,

so mu_c({0})=0. Therefore

    LH_c=mu_c-delta0,
    mu_c positive probability, mu_c({0})=0.

The source may live on multiple singular strata and need not have a classical density. Its entire mass, including any strata and infinite tails, is settled distributionally. The singular-set classification in Section6 is not needed for this source proof.

## 4. Finite-source comparison for every common radius and general L1

Here is the exact interface used for H_c. Suppose K>=0, K in L1, and

    LK=a(mu-delta0), a>0,

as a whole-space distribution, with mu a positive probability. Let M_K=integral K>0. Kernel/measure Fourier transforms use Khat(xi)=integral exp(-ix·xi)K(x)dx, so convolution introduces no extra factor. The distribution identity yields the continuous-function equality

    -psi(xi)Khat(xi)=a(muhat(xi)-1).

One may first test on Schwartz functions: K is L1 and psi grows at most linearly, giving the required absolute Fubini against a Schwartz transform. Since both sides are continuous, equality of distributions is pointwise equality. For u in H1, multiplication by uhat and the L2 Young inequalities yield

    K*Lu=a(mu*u-u) in L2.

This avoids assigning an undefined convolution to an arbitrary noncompact distribution.

For r>0, K_r(x)=r^(-7)K(x/r), and mu_r is the pushforward under z->rz. Order-one homogeneity gives

    LK_r=(a/r)(mu_r-delta0),
    K_r*Lu=(a/r)(mu_r*u-u) in L2.

The factor1/r is essential. Let f>=0 be bounded and compactly supported. Apply Section2 with kappa=lambda/M_K to get the **same** u,sigma,D for all r. On D^c, u=0 and mu_r*u>=0 a.e.; finite-measure/Lebesgue Fubini ensures this even for measures supported on low-dimensional subspaces. Therefore for each fixed r,

    K_r*f=K_r*sigma-K_r*Lu<=kappa M_K=lambda

on D^c up to a null set. Merge the null sets for positive rational r. Since dilation is L1-continuous,

    ||K_r*f-K_s*f||infinity<=||f||infinity||K_r-K_s||1 ->0,

the bound extends on the same complement to every r>0. Thus

    lambda |{sup_{r>0} K_r*f>lambda}|
      <=lambda|D|<=M_K||f||1.

For general f in L1, let f_j=min(|f|,j)1_{[-j,j]^7}, increasing to |f|. Positive convolution, Tonelli and sup_r sup_j=sup_j sup_r imply the corresponding maximal functions increase to the full one. Their strict superlevel sets increase, so measure convergence proves

    lambda |{sup_{r>0}K_r*|f|>lambda}|<=M_K||f||1.

This also gives a measurable representative of the full-scale supremum. It does not assume Linfinity continuity in r for general L1 input. “Same D” means for each bounded input and threshold across all r, not one common D for different input approximants.

Lebesgue measure is used in the mass-one convolution semigroup, self-adjoint Fourier/energy pairing, translation-invariant Lp inequalities, dx dmu-Fubini, cutoffs, volume bound kappa|D|<=integral sigma and K*sigma<=kappa integral K. These ingredients have no automatic replacement for a variable background measure. The comparison controls centred cubes with one shared radius, not seven independently selected coordinate radii.

## 5. Exact conditional g/I cost

Let Q=[-1,1]^7 and g=essinf_Q G. It is strictly positive: p_t is positive jointly continuous on compact [1,2] times Q, so its minimum there is positive and G dominates its time integral. Local L1 ensures g<infinity. This existence argument is not a useful quantitative certificate.

Assume a rigorous certificate 0<g0<=g. For 0<c<g0 set

    N=128(g0-c), K=H_c/N.

Then K>=1_Q/128 a.e., so every centred cube average with common radius r is dominated by K_r*|f|. Section3 supplies LK=(mu_c-delta0)/N, and Section4 supplies the maximal comparison. Hence C7<=integral K whenever its mass is explicitly certified.

Let q=7/6 and J=sum over the14 faces F of integral_F G(y)^q dS(y). Since G is homogeneous of degree-6, nonnegative radial integration in cube-face coordinates x=ry gives

    |{G>1}|=(1/7)J.

Indeed the Jacobian on each open face cone is r^6, and the ray condition is 0<r<G(y)^(1/6). Boundaries of face cones are null. The weak-Lq bound proves J finite before it is used in a mass formula. This argument uses the canonical extended homogeneous G and is valid with infinite values on a face-null set.

Along a ray,

    integral_0^infinity (r^(-6)G(y)-c)_+ r^6 dr
      =(6/7)c^(-1/6)G(y)^(7/6).

The formula follows by integrating to r=(G(y)/c)^(1/6), and Tonelli permits integration in y. Hence

    integral K=(6/7)c^(-1/6)J/[128(g0-c)].

Its logarithmic derivative in c equals (7c-g0)/[6c(g0-c)]. Thus the unique scalar optimum is c=g0/7. The frame law is invariant under every signed coordinate permutation: the transformation maps V to representatives in V up to column signs, and each Cauchy variable is symmetric. Therefore for

    I=integral_[0,1]^6 G(1,y)^q dy,

all14 faces have the same integral and each splits into64 equal orthants. Consequently J=14·64 I=896 I=7·128 I, and

    integral K=7^(7/6) I/g0^(7/6),
    C7<=7^(7/6) I0/g0^(7/6) if I<=I0 is certified,
    (integral K)^6=7^7 I^6/g0^7.

No point sample or floating quadrature supplies g0 or I0. The present work supplies neither full certificate and asserts no improvement of the current numerical C7 threshold. In particular, the three-/four-dimensional upper bounds are not inherited as dimension-seven upper bounds.

For essential-face reasoning, cube and boundary essential infima agree. If a positive-area face set had G<a<g, homogeneity at r sufficiently near1 would produce a positive-volume interior subset with G<g, a contradiction. Conversely G(ry)=r^(-6)G(y)>=G(y) for0<r<=1. This is an essential trace argument; lower semicontinuity alone would not extend a dense-set lower bound to all boundary points. Section6 instead supplies genuine continuity wherever G is finite.

## 6. Exact singular set: all at-most-three-direction spans

Define

    S_j=union{span(A): A subset V, |A|<=j}, j=1,2,3,

including0. These are finite unions of closed proper subspaces and have zero7-dimensional Lebesgue measure. On the face x1=1, a rank<=3 subspace has affine intersection dimension at most2, hence zero6-dimensional face measure.

### 6.1 Lower bound on the singular strata

If x in S3 is nonzero, choose a representation x=sum_{i=1}^k a_i v_i with the smallest possible cardinality k, 1<=k<=3. All a_i are nonzero. The v_i are independent: if a nontrivial relation among them existed, subtracting an appropriate scalar multiple of it from the coefficients would remove at least one nonzero coefficient, contradicting minimality. Extend these k columns to a basis B selected from V. Write

    X_t=B Z_B(t)+R_t,
    b=B^(-1)x=(a1,...,ak,0,...,0).

Choose the event A_t that every one of the57 nonbasis parameters satisfies |Z_v(t)|<=t. Its probability is (1/2)^57, independent of t. Put

    M_i=sum_{v not in B}|(B^(-1)v)_i|.

On A_t, |(B^(-1)R_t)_i|<=M_i t. For0<t<=1 and i<=k,

    f_t(a_i-(B^(-1)R_t)_i)
      >= t/[pi(1+(|a_i|+M_i)^2)].

For i>k,

    f_t(-(B^(-1)R_t)_i)>=1/[pi(1+M_i²)t].

Using the **basis density** and averaging over this positive-probability event gives the explicit positive lower bound

    p_t(x)>=c_x t^(2k-7),
    c_x=(1/2)^57/(|det B| pi^7)
       ·product_{i<=k}[1+(|a_i|+M_i)^2]^(-1)
       ·product_{i>k}(1+M_i²)^(-1)>0.

This is a lower bound on the actual spatial density, not the probability of hitting a point. The event is only used to lower-bound its convolution expectation. At0 the scaling identity gives p_t(0)=p_1(0)t^(-7)>0.

Thus the lower exponents for minimal k=0,1,2,3 are -7,-5,-3,-1. All these time integrals diverge at0, so G(x)=infinity on S3. In the k=3 case specifically,

    integral_epsilon^1 p_t(x)dt>=c_x log(1/epsilon).

This proves a logarithmic **lower divergence**, not an unproved exact asymptotic or upper logarithmic bound. Points lying in lower-rank intersections use their smallest k.

### 6.2 Rank extraction for the off-S3 upper bound

Let E be compact and disjoint from S3, d=dist(E,S3)>0. Let U={x:dist(x,E)<d/4}, so dist(U,S3)>=3d/4. Fix

    delta=d/(2m sqrt(7)).

For a parameter outcome z=(z_v), let A(z)={v:|z_v|>delta}. If rank A(z)<=3, choose at most3 independent columns spanning A(z). The sum of its large terms lies in their span, which is contained in S3. The remaining terms have Euclidean total size at most

    sum_{v not in A}|v||z_v|<=m sqrt(7)delta=d/2.

The total outcome cannot then lie in U. Consequently every parameter outcome whose spatial sum belongs to U has rank A(z)>=4. Select four **independent** large columns from A(z). This is the required independence step: having four large parameters alone is not sufficient.

There are finitely many independent4 subsets T of V. For each T, fix once and for all an extension B_T to a7-column basis (ordinary finite-dimensional basis extension applies since V spans). The event just described is covered, on U, by the finite union of events |Z_v|>delta for all v in T. We do not assume a different basis can be chosen continuously or measurably in x; there are finitely many fixed choices.

For each fixed T, the restricted one-dimensional law of a member has density

    1_{|z|>delta} f_t(z)<=t/(pi delta²).

The other3 basis densities are at most1/(pi t). Dividing by |det B_T|, and convolving all remaining57 unrestricted laws as probabilities, yields a subprobability spatial density bounded by

    t/[|det B_T| pi^7 delta^8].

Every nonzero7x7 determinant with columns from V is an integer multiple of64: subtract its first column from the other6, each now even. In particular |det B_T|>=64. Therefore, writing N4 for the number of independent4 subsets,

    p_t(x)<=N4 t/(64 pi^7 delta^8) on U,
    N4<=C(64,4)=635376.

The event-image inequality first gives the bound a.e. on the open U. Continuity of p_t promotes it to every point there. The exact finite enumeration below improves N4 to630560, but the analytic assertion does not depend on that improvement. This is a density estimate for restricted laws, not a product of event probabilities. Dependencies among all other large variables cause no loss because they are integrated as probability convolutions after an independent4 group has been selected and completed.

For t<=1 this is locally uniform O(t). For t>=1, Ct^(-7) is integrable. Dominated convergence now proves G finite and continuous on U. Since E was arbitrary, G is finite and continuous on R7 minus S3. In fact for epsilon<=1,

    0<=integral_0^epsilon p_t(x)dt
       <=N4 epsilon²/(128 pi^7 delta^8), x in E.

Together with the lower bound, this establishes the **exact** extended-value singular set {G=infinity}=S3. It replaces the false higher-dimensional extrapolation that only the64 direction lines are singular.

### 6.3 Small exact enumeration and its scope

Each three distinct V columns is independent. Otherwise their x1=1 representatives would be three collinear vertices of {-1,1}^6. A changing coordinate on the line between two cube vertices can remain in {-1,1} only at either endpoint, excluding a third vertex.

A rank3 subspace intersects the slice x1=1 in an affine2-plane. Some pair of free-coordinate projections is injective on that plane, so it contains at most4 V vertices. The analogous rank2 slice is a line and contains at most2. Thus there are64 distinct direction lines and2016 distinct rank2 subspaces. A triple-span closure contains either3 or4 V columns.

The checker processes all41664 triples. For each new span it picks two free coordinates with nonzero2x2 determinant in the difference vectors, solves for the two affine coefficients using exact integers, and tests all64 vertices against all7 coordinates. It stores the V-closure, which canonically determines the span since it contains its generating triple; all triples within a closure are marked together. It never uses a floating rank tolerance.

Actual results:

|object|count|
|---|---:|
|directions / lines|64|
|pairs / distinct rank2 spans|2016|
|triples|41664|
|rank3 spans with3 V columns|22400|
|rank3 spans with4 V columns|4816|
|all distinct rank3 spans|27216|
|dependent4 subsets|4816|
|independent4 subsets|630560|

A separate formula cross-checks the4816 count. Four vertices in an affine2-plane form a parallelogram: under an injective projection to two free coordinates, they project onto all four square corners, and every other coordinate is an affine function taking only +/-1 at these corners, hence is constant or plus/minus one projected coordinate. Opposite corners therefore have the same midpoint in all coordinates. A midpoint whose opposite pair differs in d coordinates has2^(d-1) possible unordered opposite pairs; there are C(6,d)2^(6-d) such midpoints. Thus

    dependent4=sum_{d=2}^6 C(6,d)2^(6-d) C(2^(d-1),2)=4816.

No three are dependent, so each dependent4 has rank3 and corresponds to exactly one four-vertex closure. Also41664=22400+4·4816 and27216=41664-3·4816. Only1,741,824 integer vertex-membership tests were needed; no635376-quadruple sweep was run. These counts do not prove g/I bounds and do not replace the universal basis-extension proof.

## 7. Whole-face short-time L^(7/6) budget

Define nonnegative time integrals directly,

    T_epsilon(x)=integral_0^epsilon p_t(x)dt,
    H_epsilon(x)=integral_epsilon^infinity p_t(x)dt.

Then G=H_epsilon+T_epsilon by Tonelli, including extended values. We never define T as an infinity-minus-infinity subtraction. H_epsilon is finite and continuous everywhere; T has the singular strata of Section6 but a finite face Lq norm.

### 7.1 One-large-variable spatial-density estimate

Every first coordinate v1=1. If the spatial sum satisfies x1>1-eta, some one of the64 coefficients has |Z_v|>=(1-eta)/64. For each v choose the previously verified determinant512 basis containing v. On its large-parameter event,

    f_t(z)1_{|z|>=(1-eta)/64}
      <=64² t/[pi(1-eta)²].

The six complementary basis densities are <=1/(pi t). The Jacobian is1/512. Convolving the other57 variables as probabilities gives a subprobability spatial-density upper bound64²t^(-5)/[512pi^7(1-eta)²]. Sum all64 possible large choices; event overlaps are harmless in this upper bound. The true density therefore obeys, a.e. in this open halfspace and hence everywhere there by continuity,

    p_t(x)<=512 t^(-5)/[pi^7(1-eta)²].

Evaluate at x1=1 and let eta down0 to get

    sup_y p_t(1,y)<=C_f t^(-5), C_f=512/pi^7.

This argument explicitly works with the density of a restricted law. A fixed face is a zero-probability slice; multiplying event probabilities would not prove its density bound. The direction columns in the determinant are unnormalised and have weight1; scale t is paid once in each Cauchy density.

### 7.2 The same density's exact normal marginal

The first coordinate of X_t is a sum of64 independent scale-t Cauchy variables, so its normal marginal has scale64t. To justify the fixed-slice identity without relying on an arbitrary a.e. representative, first integrate the explicit basis product q_B,t(a,y) over y. Changing variables and using that the first row of B is all1 gives the7-variable Cauchy convolution density7t/[pi(a²+49t²)] for every real a. Convolution with the other57 variables and nonnegative Tonelli then gives, for every a,

    integral_R6 p_t(a,y)dy=64t/[pi(a²+4096t²)].

In particular at a=1 this is <=64t/pi. This is a marginal of the same p_t used in the sup estimate, not a product of separate-coordinate processes.

### 7.3 Interpolation, time integration and orthant factor

Set q=7/6. The elementary inequality p^q<=||p||infinity^(q-1)p, together with the two preceding bounds, gives

    ||p_t(1,·)||_(Lq(R6))
      <=C_f^(1/7)(64/pi)^(6/7)t^(1/7)
         (1+4096t²)^(-6/7)
      <=C_f^(1/7)(64/pi)^(6/7)t^(1/7).

The exponent is exactly [(-5)(q-1)+1]/q=1/7. Apply Minkowski first on finite time/space truncations and then positive monotone limits:

    ||T_epsilon(1,·)||_(Lq(R6))
       <=(7/8)C_f^(1/7)(64/pi)^(6/7)epsilon^(8/7).

The law is invariant under reflection in each of the six free coordinates, as is T_epsilon. Hence the qth integral on the positive orthant is2^(-6) of the whole-plane integral; the **norm** factor is

    (2^(-6))^(1/q)=2^(-36/7).

Since [0,1]^6 is contained in that orthant,

    ||T_epsilon(1,·)||_(Lq([0,1]^6))<=E_epsilon,
    E_epsilon=2^(-36/7)(7/8)(512/pi^7)^(1/7)
                (64/pi)^(6/7)epsilon^(8/7)
             =(7/8)2^(9/7)pi^(-13/7)epsilon^(8/7).

This is valid for every epsilon>0, covers the entire face including all three-rank singular sets and all spatial tails in the larger plane, and deletes no caps. It does not assert that G itself is bounded near S3.

If a continuous certificate proves integral_[0,1]^6 H_epsilon(1,y)^qdy<=B, then

    I<= (B^(6/7)+E_epsilon)^(7/6).

The budget expression E_epsilon has rigorously enclosed values in `results.json`; its lower endpoint is not a lower bound on the unknown actual short-time norm. Only an upper endpoint can be used for a certified total I bound. The fixed expressions were evaluated with256-bit Arb, with no Green-node integration.

## 8. Reproduction, provenance and remaining proof obligations

From this independent clone run:

    /Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python -B cube-dim-7-bounds/upper/round001/theory/check.py

The finite parameters and stopping conditions were recorded before that run in `registration.md`. The first actual run completed with exit0, no failed assertion and no repair; elapsed script time was2.244968750019325 seconds. `run.stdout` records its displayed results; `run.stderr` is empty. Source hashes were taken before and after and all matched. Python, SymPy, python-flint and precision versions are recorded in `results.json`. SymPy/Fraction-style exact integer arithmetic and Arb are shared libraries, so this is not a claim of independent arithmetic kernels.

All writes were restricted to this owner's directory in the independent clone. No commits, pushes, MATLAB calls, recursive agents, prior node runs or optimization were made. The read-only evidence is not copied into this deliverable; `sources.json` provides paths and SHA256 provenance. `run_receipt.json` records this execution, and `manifest.json` whitelists the final files.

The remaining substantive obligations for a useful C7 upper bound are a globally valid useful g0 and a certified full I0 for this exact common-time frame, plus a strict exact comparison with the current numerical threshold. The analytic source/comparison proof already handles general L1 and all common radii and does not need a future pointwise Green formula. The S3 classification permits continuous extension only off S3; a future candidate formula must settle every artificial denominator boundary there, and must not quietly remove neighbourhoods of singular planes. Alternatively, H_epsilon can be certified continuously and this whole-face Lq budget can pay all omitted time. This round stops at the registered interface and finite structural checks.
