# R071 · A finite second-profile Cantelli cut on the genuine cube history

Status: full proposed finite theorem, pending independent review. Every
claim in this note uses the one fixed tuple below, without a parameter
search, a Gaussian or Berry–Esseen approximation, or a change of the
original maximal operator. A local or scalar check alone is not a C_n
theorem. No frozen earlier source or another agent's file is changed.

## 1. Fixed data and proposed assertion

Let L=log n>=L0=8192, x=91L/100, t=n/x, beta=1+1/t and

    s=3/2, qmin=2/3, qmax=91/100,
    delta=2^-19, rlo=1+delta, rhi=273/200-delta,
    h=1/(256L), zeta=n^-3, epsilon=zeta t,
    kappa=beta-zeta, sigma=9/100, chi=9/200,
    theta=exp(-chi x), k=nh theta/rhi.

The tolerances are

    eB=eN=eT=2^-19,
    eJ=eRad=eS=eD=2^-20,
    epsA=1/65536, epsRet=14/25, Ret=11/25.

Define exactly

    P0=1-eB-6/256,
    Ulo=P0(1-eN)/(1+eJ), Uhi=(1+eN)/[P0(1-eJ)],
    M=2(1+eRad)Uhi,
    b=2 Ret (1-eRad)Ulo-eS-eT,
    mL=311/1000, Gamma=mL b/M, gamma=13/100.

The proposed conclusion is the complete original-operator bound

    C_n<=F(A_n(t),B_n(t)-gamma n)
       <n log n-n log log n-(7/200)n,                   (1)
    F(A,B)=(sqrt A+sqrt B)^2.

This is for all integer n with log n>=8192, arbitrary L1 input, every
level and every positive common witness radius. A_n,B_n are the original
cover and corner-envelope quantities, not a substituted model cost.

The exact single-tuple verifier proves 43 finite rational gates. In
particular (1-eD)Gamma>gamma, and

    (1-epsA)log(rhi/rlo)>mL,
    -log(1-sigma)<189/2000,
    gamma-189/2000-7/200=1/2000.                        (2)

## 2. Full bands, all active guards and the final bounded kernel

Use the original corner envelope

    D={||z||_infinity<t: some z_i<1-beta m(z)},
    m(z)=max(1,||z||_infinity),
    d lambda=2^-n m(z)^-n 1_D dz.

Write lambda0 for the same density on full annuli before D restriction.
At a shell radius rho>1 each signed maximal face has measure
d rho/(2rho), with n−1 independent Uniform[-rho,rho] free coordinates.

Let A0={w:r=||w||_infinity in [rlo,rhi]}, let Pfull be the full source
band q=||z||_infinity/t in [qmin,qmax], and put P=D intersect Pfull.
Matching is |r-sq|<h. Throughout the FULL source band define

    rho=tq, R=sq-h, a0=(1-epsilon)q,
    g=R-a0=(s-1+epsilon)q-h,
    ell_R(z_i)=min(1,(R+kappa rho-|z_i|)/(2R)),
    G(z)=product_i ell_R(z_i), bsoft(z)=min(1,theta/G(z)).

These definitions are made before removing the source cap complement.
They have R>0,a0>0,G>0 and 0<bsoft<=1. Extend bsoft=0 off Pfull.
Set arad(w)=r/rhi on A0 and zero off A0. After A is defined below,

    K(w,z)=1_A(w)1_P(z)1_{|r-sq|<h}
      1_{||z-w||_infinity<kappa rho} arad(w)bsoft(z).    (3)

Write R_lambda(w)=integral K(w,z)d lambda(z) and
C_lambda(z)=integral K(w,z)d lambda(w).

For target r the matching q endpoints and exact baseline are

    q_-=(r-h)/s, q_+=(r+h)/s, qb=q_+.

The entire interval lies in [qmin,qmax], since 2h<delta. At qb,
R=r; on matching R lies in [r−2h,r]. Active sources satisfy
R>rlo−2h>1. Unmatched edge sources can have R<1 but have zero final
column; no R>=1 local lemma is applied there. The finite gates give

    q>=2/3, 1<R<=r<3/2, rho>=t qmin>4rhi>4,
    h<1/1000, epsilon qmax<=2h/3, epsilon<1/2,
    0<g<=R/3, kappa>1, beta<2.                         (4)

Indeed g−R/3=epsilon q−2h/3<=0. The positive zeta rho term is
retained. Both annuli lie within norm<t, and A0 and P are disjoint.

## 3. Coordinate profiles and explicit baseline moments

At fixed target coordinate |w|<=r define, all before source-cap removal,

    J(w;q)=(2rho)^-1 integral 1_capture/ell_R dz,
    J1(w;q)=(2rho)^-1 integral 1_capture phi/ell_R dz,
    J2(w;q)=(2rho)^-1 integral 1_capture phi^2/ell_R dz,
    capture={z in [-rho,rho]: |z-w|<kappa rho},
    phi(z)=-log ell_R(z).

Thus J2 is a second integrand profile, NOT the square of J1. Put
Q(R,u)=−log(1−u/(2R)), m=(|w|−a0)_+, v=(g−m)_+, and
F(R,u)=−2R log(1−u/(2R))−u. Direct integration gives

    J=1+[F(R,g)+F(R,v)−m]/(2rho),
    J1=R[Q(R,g)^2+Q(R,v)^2]/(2rho),
    J2=R[Q(R,g)^3+Q(R,v)^3]/(3rho).                   (5)

Central deleted segments when m>g have phi=0, so the two profile
formulas hold also when |w|>R. Positive-part onsets are handled by
local Lipschitz continuity and absolute continuity; pointwise classical
differentiability at every onset is not asserted.

The relative-coordinate bounds on (4) are

    |partial_q log J|<3/rho, |J−1|<1/(3rho),
    0<=J1<=1/rho, |partial_q J1|<1/(2rho),
    0<=J2<=1/rho, |partial_q J2|<=1/rho.               (6)

The first two follow from the accepted relative lemma (using its actual
3/rho bound, not the older loose 256/rho). The first-profile derivative
is the stronger bound proved in the accepted R066 local lemma. For the
second profile the new tail formula (5) proves an even stronger bound:
Q<=1/5, |Q_R|<=1/(5R), Q_u<=3/(5R), and, writing
S3=Q(R,g)^3+Q(R,v)^3,

    S3<=2/125, |S3'|<=63/(250R),
    (R/(3rho))'=(h/q)/(3rho),
    |J2'|<1/(10rho).

The weaker 1/rho in (6) is deliberately used in the finite guards.
All these derivative statements are a.e./AC on the full matching
interval, and are unconditional in the target coordinate value.

At qb, R=r, let W be a full TARGET Uniform[-r,r] coordinate. Then
E_W J=1, and

    -6g^3/rho^2<=E log J<=0,
    Var(log J)<=12g^3/rho^2.                            (7)

The pinned target log term lies in [−g/rho,0]. Fubini cancels the
inverse ell from the two new profile means. With a=1−g/(2R)>=5/6,

    E_W J1=E_Z phi=(2R/rho)[1−a+a log a],
    E_W J2=E_Z phi^2
           =(2R/rho)[2−a((log a)^2−2 log a+2)].        (8)

Here Z is FULLSOURCE Uniform[-rho,rho]; this identity does not identify
the target and source inverse laws. Define their limiting coefficients

    mu1=1/2−(5/2)log(6/5)<221/5000,
    c2=3[2−(5/6)(Q0^2+2Q0+2)], Q0=log(6/5).

The rational log certificate encloses c2 strictly between .0052893
and .0052894; decimal endpoints here mean the exact corresponding
terminating rational numbers. For the mean of the whole TARGET-face
J2 sum, let n>=N and use baseline

    a−5/6=h/(3R)−epsilon q/(2R) in [0,h/(3R)].

Since the derivative of integral_a^1 log^2 u du is −log^2 a>=−1/25,
and R/q=s−h/q, its free-coordinate mean coefficient is between
c2−(c2+1/25)h and c2. Omitting one free coordinate and adding the
nonnegative deterministic pin gives safe whole-face bounds

    (1−1/N)[c2−(c2+1/25)h]
       <=E sum J2/x<=c2+1/(N qmin).                   (9)

In particular these coefficients lie in [.005275,.005295] at and
above the base gates. The weak pointwise 1/rho in (6) yields
Var(sum J2)<=n/rho_b^2; sharper small constants are unnecessary.

## 4. Six target conditions and finite good mass

For fixed target w define

    N_i(q)=1+1_{|w_i|<a0(q)}, Nsum=sum_i N_i,
    E_±=n+(n−1)a0(q_±)/r, ellpin=(1+a0/R)/2.

On a full signed target face the pin has N_i=1 and the other n−1
indicators are independent Bernoulli with variance sum<=n/4. Let A
be A0 intersected with these SIX events:

1. |sum_i log J(w_i;qb)|<=eB.
2. |Nsum(q_-)−E_-|<=eN n/4.
3. |Nsum(q_+)−E_+|<=eN n/4.
4. w belongs to D.
5. sum_i J1(w_i;qb)<=(221/5000+1/100000)x.
6. .00527x<=sum_i J2(w_i;qb)<=.0053x.

Every stated decimal is an exact rational, and each event concerns the
ambient full-target-face law. They are not imposed as a distribution
assumption on actual increments, and independence between these six
events is not required.

For explicit bounds set N(L)=L^64/64!<=e^L=n. Uniformly on each
signed target face the log baseline mean and variance are bounded by

    ML=6(1/8)L^2/(qmin^2 N)+(1/2)L/(qmin N),
    VL=12(1/8)L^2/(qmin^2 N).

The finite gates ML<=eB/2 and 4VL/eB^2<=epsA/8 prove the first bad
probability<=epsA/8. Each count failure is at most 4/(eN² N), also
at most epsA/8. The positive pinned target cap-complement fraction
is at most (3/4)^(n−1); negative pinned faces contribute zero. The
verified target-fraction gate and n>=75 make it at most epsA/8.

For condition 5 the accepted profile estimate gives failure at most
324000000/n once n>=9000. This follows from (8), the deterministic
target pin bound and independent free target coordinates; it is at
most epsA/8 using N. For condition 6, (9) puts the mean at distance
at least 1/200000 times x from both excluded endpoints. The variance
bound after (9) and n/rho_b=x/q_b give failure at most

    1/[N qmin^2(1/200000)^2]<=epsA/8.

The six bad probabilities sum to at most 6epsA/8<epsA. Therefore

    lambda(A)>=(1−epsA)n log(rhi/rlo)>mL n.             (10)

All finite inequalities used here occur in the independent exact
43-gate author certificate. Section 7 proves their extension beyond
the base point, so this is not a finite-sampling argument.

## 5. Whole-matching profiles and per-source-face Cantelli

The count calibration has the correct multiplicative sign:
ellpin'=−(1−epsilon)h/(2R²)<0 and Nsum is nondecreasing. Thus the
ratio Nsum/(2n ellpin) is nondecreasing. Its endpoint mean shortfall
at q− is at most 2h+1/n and its mean upper ratio at q+ is at most
one. A count fluctuation changes the ratio by at most eN/4. The
finite gate 2h+1/n<=3eN/4 therefore gives

    1−eN<=Nsum/(2n ellpin)<=1+eN

over the entire matching interval. The strong first derivative in
(6), n/rho<=3L/2 and interval width4h/3 give

    |sum log J(q)−sum log J(qb)|<=6Lh=6/256.

Hence P0<=product J_i<=1/P0. Also the individual gate
x/(3qmin N)<=eJ gives every J_i in [1−eJ,1+eJ] pointwise,
without a union bound. The complete inverse source-shell coefficient is

    U(w,rho)=[product_i J_i/(2ellpin)]sum_i N_i/J_i,
    Ulo n<=U(w,rho)<=Uhi n.                            (11)

The remaining derivatives in (6) give deterministic transport

    |sum J1(q)−sum J1(qb)|<=xh<=x/100000,
    |sum J2(q)−sum J2(qb)|<=2xh<=x/1000000.

Therefore for every w in A and EVERY matched source radius,

    sum J1<=2211x/50000,
    .005269x<=sum J2<=.005301x.                        (12)

Now fix one CAPTURED SIGNED SOURCE FACE j in its full inverse-capture
law, before corner-cap removal. Its pin contributes Qpin<=log(6/5)
deterministically to V=−log G; all other coordinates are independent,
with means J1_i/J_i and second moments J2_i/J_i. Consequently

    E(V|face j)<=.183+[2211/(50000(1−eJ))]x,
    Var(V|face j)<=sum_(i!=j)J2_i/J_i
                   <=.005301x/(1−eJ)<.00531x.           (13)

For the upper bound only, the omitted source pin can be added back
because J2 is nonnegative. The lower half of the two-sided profile
(12) is true but NOT USED in (13): Cantelli here requires only a
variance upper bound. No lower variance, subtraction of mean squares,
or unproved mixture variance is used.

At L>=8192 the exact finite gates give

    chi−2211/[50000(1−eJ)]−.183/x>1/1330,
    (.00531)1330^2/x<=14/11.

Apply Cantelli separately within EACH face. If its variance is v and
its positive threshold gap is u, then

    P(V>chi x|face j)<=v/(v+u²)<=14/25.

Only afterwards mix these uniformly bounded probabilities over
source faces. This retains at least Ret=11/25 of the complete inverse
row, without claiming that distinct face means coincide. Since bsoft
divided by theta/G equals one on G>=theta and lies in [0,1] elsewhere,
the full SOFT shell coefficients are between Ret theta Ulo n and
theta Uhi n.

## 6. Complete calibrated degrees, source cap and original tau

The exact source logarithmic width is log((r+h)/(r−h)). The radial
factor r/rhi and 1<=atanh(u)/u<=1/(1−u²), u=h/r, give the bounds

    (1−eRad)2h/rhi <=(r/rhi)log((r+h)/(r−h))
                                  <=(1+eRad)2h/rhi.

The full-source raw row is thus between 2Ret(1−eRad)Ulo k and Mk.
For a fixed source and target radius r>=R, the captured target sign
count T_i obeys T_i/(2ell_i(r))<=1: if ell_i<1 there is one captured
sign, while if ell_i=1 there are at most two. Summing ALL target faces
gives coefficient<=n product ell_i(r)<=nG(z). The equality
(dr/r)arad=dr/rhi and matched target width<=2h prove the complete
column bound C_lambda<=2nh theta/rhi=2k. Inactive columns are zero.
Thus M>2 controls both ambient degrees, after any restrictions.

The full source band mass is at most Csrc n, where
Csrc=(qmax−qmin)/qmin. Its original cap-complement mass is at most

    Csrc n exp[−(1−1/n)eta x], eta=(1−qmax)/(2qmax)=9/182.

Because both kernel weights are <=1, this UNWEIGHTED mass pays its
removal from each raw row. The exact rate relations are

    (eta−chi)x=(81/20000)L,
    [sigma/2−chi(1−sigma)]L=(81/20000)L.

The finite-n correction exp(eta x/n) is at most 2, by eta x/N<=1/2.
Hence cap loss divided by k is at most

    2Csrc rhi 256L exp[−(81/20000)L]<=eS.              (14)

The original scalar cap identity is B_n=core+n log t−(n/2)Icap,
with core>=0 and 0<=Icap<=log t. It gives B_n>=(n/2)log t.
For L>=8192, exp(L/2)>=(L/2)^2/2>=L>=x, so log t>=L/2.
Thus B_n>=nL/4>gamma n. Also B_n<=1+n log t<=nL, and
A_n=beta^n<=e^x=n exp(−sigma L). Choose the original ONE threshold

    tau=sqrt[A_n(B_n−gamma n)].

It is positive and obeys

    tau/k <=rhi 256L^(3/2)exp[−(81/20000)L]<=eT.        (15)

Equations (14)–(15) use the actual normalization k with rhi; no factor
is discarded. Final row minus the SAME global tau is at least b k,
and both final degrees are at most Mk.

## 7. The finite 8192 gates and the all-tail proof

The author verifier uses Fraction arithmetic for N(L)=L^64/64! and
the FULL prefix S64(y)=sum_(j=0)^64 y^j/j!, so exp(y)>=S64(y).
It does not use a single Taylor term for the delicate cap/tau bounds.
At L0=8192, sqrt L0<91, and it verifies

    2Csrc rhi 256L0/S64((81/20000)L0)<=eS,
    rhi 256L0*91/S64((81/20000)L0)<=eT.                 (16)

These bounds extend to every L>=L0. For y>=2,
y S64'(y)−S64(y)>=−1+y²/2>=0. For y>=4,
y S64'(y)−(3/2)S64(y)>=−3/2−y/2+y²/4>=0.
Thus L/S64(cL) and L^(3/2)/S64(cL) decrease on the required tail.
The base y is greater than4. This proves (14)–(15) for the whole tail,
including the irrational square root; its rational upper91 is used
only to check the true base expression before monotonicity.

All other geometry/probability upper bounds are constants, negative
powers of L or positive sums thereof, after N=L^64/64! is inserted.
In particular epsilon qmax/h=O(N^-2), ML and VL decrease, individual
J errors and count/profiles bad probabilities decrease, and the count
and radial errors decrease. The J2 mean lower envelope in (9)
increases and its upper envelope decreases. The positive Cantelli gap
rate increases (only .183/x varies) and its ratio upper decreases.
The lower source-radius bounds N/x increase. The product loss6Lh
is exactly constant; h and trim guards improve. Therefore checking
these explicit bounds at8192 covers every larger L. No numerical
sampling over a dimension tail is used.

## 8. Same actual finite history and true missing mass

Use the original largest-eligible-radius history with real disjoint
increments I_i. For algorithmic dilation mathfrak a>=8n^4 normalize
at any positive load-support point by its last minimum positive-radius
contributor and reflect its center into the positive orthant. Write
alpha=mathfrak a/(mathfrak a+1), ccap=(mathfrak a beta−1)/(mathfrak a+1).
The actual uncaptured-increment measure restricted to D is dominated
by lambda_a=alpha^-n lambda. The full finite-cap deficit is no smaller
than this restricted deficit, since D lies in the finite cap by

    alpha−ccap m−(1−beta m)
       =[(beta+1)m−1]/(mathfrak a+1)>0.

For a target actual item j meeting A, its uncaptured-center inequality
gives ccap r_j<=||u−w||<alpha+rhi<3; since ccap>1, r_j<3.
For a source point z=y/(mathfrak a+1), the common support and y in Q_i
give ||z||<r_i. Hence r_i>=rho>4, and source i precedes the SAME
target j. The displacement guard

    alpha beta−kappa−2/(mathfrak a+1)
       =zeta−(beta+2)/(mathfrak a+1)>0

holds for all mathfrak a>=8n^4. Included pairs really capture that
target center in its own eligibility test. Integrating only included
increment pieces costs at most their original item weights. Both arad
and bsoft are <=1, so actual weighted rows are <=arad(w)tau<=tau.
No ambient law is imposed on actual mu_a.

Let E be the limiting ambient energy of K, ca=alpha^-n>=1, and
nu_a=lambda_a−mu_a. The energy difference identity, disjoint A,P,
and scaled degree bound ca Mk give

    Delta_a >=[ca² E−tau ca lambda(A)]/[ca Mk]
            >=[E−tau lambda(A)]/[Mk]
            >=lambda(A)b/M>=Gamma n>gamma n.          (17)

This is a missing-mass statement for every original finite history,
uniformly in algorithmic a. The ca²/ca factors are retained; no limit
of an actual selector or history is presumed.

## 9. General L1, every positive common radius and complete root cost

The original finite cover/load/witness argument, with captured budget
tau and final self weight<=mathfrak a^-n, now yields

    C_n <=[1+(beta+1/mathfrak a)^n/tau]
       [tau+mathfrak a^-n+alpha^-n J_mathfrak a(t)−gamma n].

For arbitrary L1 input use its absolute value. Every point of a strict
maximal level set has a true positive witness radius; its open small
cube covers that point, and finitely many cover every compact subset.
The selected dilates stay true witnesses. The above normalization works
for every positive terminal radius. Inner regularity gives the full
level set; no compact-support or scale-truncation hypothesis is added.
For fixed n let mathfrak a tend to infinity only in the uniform scalar
bounds, obtaining C_n<=F(A_n,B_n−gamma n).

The original scalar estimates give, with T=n log n−n log log n,

    B_n<=T−n log(1−sigma)+1,
    A_n+2sqrt(A_n B_n)<=n e^-sigma L+2n sqrt L e^-sigma L/2.

The exact price bound in (2) leaves slack1/2000 after paying7/200.
At8192 the verifier checks

    1/N(L0)+1/S64(sigma L0)+182/S64(sigma L0/2)<1/2000.

This bounds the true scalar remainder
1/n+e^-sigma L+2sqrt L e^-sigma L/2 at the base; the same prefix
monotonicity proves its decrease thereafter. Hence the complete root
cost is strictly less than T−7n/200. This proves the proposed (1),
subject to independent acceptance of the analytic and gate chain.

## 10. Evidence and dependency scope

`verify.py` was actually run once for this new fixed tuple and returned
PASS for43 exact rational gates. `verification.json` records the
results. It uses positive64-term atanh logarithm enclosures with a
geometric tail and exact full exponential prefixes. These computations
are not production searches or independent evidence of a geometric
law; all law/measure transitions are specified above.

Imported local dependencies are the accepted R064 relative-coordinate
lemma c60d48d116863c9c9cac177e19a9fa4158f12fd8dcef0d5239a5ad0a5487d9f1
for (6)–(7), and the accepted R066 first-profile lemma
67bdf01d9aa65c70bf058df08b527b904e954d19f3327ecbb43e5b5c316b5030
for the strong J1 derivative and explicit target failure probability.
Their abstract active guards, not their old h/tuple threshold, are
checked anew here. The new R070 second-profile lemma is being
independently derived; (5), (8)–(9) and the safe source variance upper
are displayed here so its identity can be bound without hiding a
lower-variance or mixture assumption. This source does not require a
source variance lower bound or any normal approximation.
