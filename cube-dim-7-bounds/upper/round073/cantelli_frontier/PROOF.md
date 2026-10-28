# R073 candidate — a newly certified finite Cantelli frontier tuple

Status: **complete proposed finite theorem, pending independent whole-chain review**. A numerical frontier model is not activated as a C_n theorem. This note registers one exact tuple, recalculates its finite guards, and proves their all-tail extensions and operator transitions. No R071 fixed parameter or threshold is imported as a verified gate.

The proposed conclusion, for every integer n with L=log n>=8192, is

    C_n<=F(A_n(t),B_n(t)-(291/1000)n)
       <n log n-n log log n-(47/250)n,                   (1)
    F(A,B)=(sqrt A+sqrt B)^2.

The operator is the original centered axis-parallel common-radius cube maximal operator on Lebesgue measure, for arbitrary L1 input, every level and all positive witness radii. No claim below covers log n<8192.

## 1. Fixed rational candidate and new exact certificate

The original candidate fields were read as rationals from `frontier4_archive/frontier4_recheck.json`, not reconstructed from decimal outputs. Its proposed L=40448/5 is replaced by the NEW base L0=8192. Fix

    s=3/2, qmin=2/3, qmax=11569/12800,
    sigma=25/256, chi=63/1280,
    x=(1-sigma)L=(231/256)L, t=n/x,
    beta=1+1/t, zeta=n^-3, epsilon=zeta t, kappa=beta-zeta,
    delta=2^-19, rlo=1+delta, rhi=s qmax-delta,
    h=1/(768L), theta=exp(-chi x), k=nh theta/rhi.

In particular qmax is NOT 1-sigma. The source-cap and tau rates in Section 7 are different. The errors and retained fraction are

    eB=eN=eT=2^-19,
    eJ=eRad=eS=eD=2^-20,
    epsA=1/65536, Ret=9717/10000, epsRet=283/10000.

Define exactly

    P0=1-eB-6/768,
    Ulo=P0(1-eN)/(1+eJ), Uhi=(1+eN)/[P0(1-eJ)],
    M=2(1+eRad)Uhi,
    b=2 Ret(1-eRad)Ulo-eS-eT,
    mL=3043/10000, Gamma=mL b/M,
    gamma=291/1000, price=257/2500, c=47/250.

The actual new one-process Fraction run verifies 53 exact gates, including

    0<P0<1, M>2, b>0,
    (1-epsA)log(rhi/rlo)>mL,
    0<gamma<(1-eD)Gamma<Gamma<1,                         (2)
    -log(1-sigma)<price,
    gamma-price-c=1/5000.                              (3)

The complete budget (1-eD)Gamma is about .291082178; this decimal is only a display of the exact rational in the receipt. The chosen gamma and c are the fixed fractions above. Registration and its coefficient-only amendment preceded the complete finite-gate run and remain preserved.

## 2. Full bands, bounded kernel, and new active guards

Use the original corner envelope

    D={||z||_infinity<t: some z_i<1-beta m(z)},
    m(z)=max(1,||z||_infinity),
    d lambda=2^-n m(z)^-n 1_D dz.

Before cap restriction write lambda0 for the same density on full annuli. Each signed maximal face at radius rho>1 has weight d rho/(2rho) and n-1 independent Uniform[-rho,rho] free coordinates.

Set A0={w:r=||w||_infinity in [rlo,rhi]}. Let Pfull be the FULL source band q=||z||_infinity/t in [qmin,qmax], and P=D intersect Pfull. Matching is |r-sq|<h. On ALL of Pfull define

    rho=tq, R=sq-h, a0=(1-epsilon)q,
    g=R-a0=(s-1+epsilon)q-h,
    ell_R(z_i)=min(1,(R+kappa rho-|z_i|)/(2R)),
    G(z)=product_i ell_R(z_i), bsoft(z)=min(1,theta/G(z)).

This includes sources outside D and unmatched edge sources. One has R>0,a0>0,G>0 and 0<bsoft<=1 on the finite domain established below. Extend bsoft by zero off Pfull. Set arad(w)=r/rhi on A0 and zero off A0. After A is specified, the FINAL kernel is exactly

    K(w,z)=1_A(w)1_P(z)1_{|r-sq|<h}
       1_{||z-w||_infinity<kappa rho} arad(w)bsoft(z).    (4)

Write R_lambda(w)=integral K(w,z)d lambda(z) and C_lambda(z)=integral K(w,z)d lambda(w). All fields and sets are Borel.

Put N(L)=L^64/64!<=e^L=n. The new base gates, extended in Section 8, give for ALL L>=8192

    2h<delta, h<=2^-20<1/1000,
    qmax/(xN^2)<=h/s, 1/(xN^2)<1/2,
    Nqmin/x>4rhi, N/x>rhi,
    x/N<1, xN^2>1, 2<=x<=L.                            (5)

They imply epsilon qmax<=h/s, beta<2, kappa>1 and rho>=tqmin>4rhi. The COMPLETE matching interval and baseline are

    q_-=(r-h)/s, q_+=(r+h)/s, qb=q_+.

They lie in [qmin,qmax], and R(qb)=r. An ACTIVE source matches some r in [rlo,rhi], and then

    R>rlo-2h>1, R in [r-2h,r],
    g-R/3=epsilon q-2h/3<=0, 0<g<=R/3.                 (6)

Inactive sources can have R<1; their final columns are zero and no R>=1 local lemma is used there. Both annuli lie within norm<t, and A intersect P is empty: targets have r<3/2 and sources rho>4. All active guards are therefore checked anew for h=1/(768L); the different h_c interval in the older LOCAL finite-gate section is not invoked.

## 3. Three exact profiles and their probability laws

For fixed target |w|<=r define, before source-cap removal,

    J_j(w;q)=(2rho)^-1 integral_capture phi(z)^j/ell_R(z) dz,
    j=0,1,2, J_0=J,
    capture={z in [-rho,rho]: |z-w|<kappa rho},
    phi=-log ell_R.

J_2 is a second INTEGRAND profile, not J_1 squared or the second moment of the whole source-face sum. With Q(R,u)=-log(1-u/(2R)), m=(|w|-a0)_+, v=(g-m)_+, and F(R,u)=-2Rlog(1-u/(2R))-u, tail integration gives

    J=1+[F(R,g)+F(R,v)-m]/(2rho),
    J_1=R[Q(R,g)^2+Q(R,v)^2]/(2rho),
    J_2=R[Q(R,g)^3+Q(R,v)^3]/(3rho).                  (7)

Extra central deletion when m>g has phi=0, so these formulas cover |w|>R. Positive-part onsets can have genuine kinks; all profiles are locally Lipschitz and absolutely continuous on the closed matching interval.

The frozen local proofs cited in Section 11, on the actual active guards (5)--(6), supply the strong a.e./AC bounds

    |partial_q log J|<3/rho, |J-1|<1/(3rho),
    0<=J_1<3/(50rho), |partial_q J_1|<1/(2rho),
    0<=J_2<1/(125rho), |partial_q J_2|<1/(10rho).       (8)

Here R'=3/2, g'=1/2+epsilon and a0'=1-epsilon. For example Q<=1/5, |Q_R|<=1/(5R), Q_u<=3/(5R); differentiating BOTH Q^3 terms AND R/(3rho) proves the last derivative bound. These are relative-gap bounds, not an extension obtained by relabeling an absolute g<=1/3 domain.

At qb, R=r. Under a full TARGET W~Uniform[-r,r], Fubini gives E J=1 and

    -6g^3/rho_b^2<=E log J<=0,
    Var(log J)<=12g^3/rho_b^2.                          (9)

The pinned target log term is deterministic in [-g/rho_b,0]. This is distinct from the captured inverse SOURCE law below.

Let Q0=log(6/5). Fubini for the two profiles, canceling one inverse ell by the target capture fraction, gives

    E_W J_1=E_fullsource phi=(2R/rho)[1-a+a log a],
    E_W J_2=E_fullsource phi^2=(2R/rho)K_2(a),
    a=1-g/(2R)>=5/6,
    K_2(a)=2-a[(-log a)^2+2(-log a)+2].

Their limiting coefficients are

    mu1=1/2-(5/2)Q0<221/5000,
    mu2=3K_2(5/6)=1-5Q0-(5/2)Q0^2.                   (10)

The new certificate encloses Q0,mu1,mu2 by exact rational logs. The upper direction for mu1 uses a LOWER bound on Q0. The exact multiplicative structure further gives

    a-5/6=h/(3R)-epsilon q/(2R) in [0,h/(3R)],
    R/q=s-h/q.

Since K_2'(a)=-(log a)^2 and Q0<1/5, the full baseline target-face sum S_2=sum_i J_2(w_i;qb) satisfies

    E S_2/x>=(1-1/n)[mu2-h(mu2+1/25)],
    E S_2/x<=mu2+1/(250nqmin).                         (11)

The final term pays the deterministic target pin, less than 1/(250rho_b). The lower mean uses the exact finite-h structure, not only the relative upper-gap guard.

## 4. Six target conditions and new finite probabilities

On each full target face, n-1 coordinates are independent Uniform[-r,r], and one is pinned. Define

    N_i(q)=1+1_{|w_i|<a0(q)}, Nsum=sum_i N_i,
    E_+/-=n+(n-1)a0(q_+/-)/r.

Let A be A0 intersected with these SIX events:

1. |sum log J(w_i;qb)|<=eB.
2. |Nsum(q_-)-E_-|<=eN n/4.
3. |Nsum(q_+)-E_+|<=eN n/4.
4. w belongs to D.
5. sum J_1(w_i;qb)<=(221/5000+1/100000)x.
6. (527/100000)x<=sum J_2(w_i;qb)<=(53/10000)x.

No independence between these events is assumed. Their laws are full-target ambient laws, not a model for actual increments.

As g<1/2, (9) bounds the log-sum mean and variance by

    ML=(3/4)L^2/(qmin^2 N)+(1/2)L/(qmin N),
    VL=(3/2)L^2/(qmin^2 N).

New gates ML<=eB/2 and 4VL/eB^2<=epsA/8 pay event 1. Since a0(q)<r, the pinned target contributes N_i=1. The other N_i-1 are independent Bernoulli(a0/r), so each count variance is <=n/4 and each failure probability is <=4/(eN^2N)<=epsA/8.

For event 4, negative pinned faces are excluded from the original cap complement. On positive faces its free fraction is 1+1/(2t)-1/(2r). The new gates

    1+x/(2N)-1/(2rhi)<=3/4,
    N>=75, (3/4)^74<=epsA/8

pay this cap condition uniformly in r.

For event 5, (10) gives free-coordinate mean <=(221/5000)x. The deterministic target pin J_1<3/(100rho_b) is <=x/(2*100000) at n>=9000. The bound J_1<3/(50rho_b) gives variance <=9n/(2500rho_b^2); Chebyshev at the remaining x/(2*100000) gap yields failure

    <=324000000/n<=324000000/N<=epsA/8.                (12)

For event 6, the NEW rational log bounds in (11) and actual h give

    (528/100000)x<E S_2<(5295/1000000)x.               (13)

The base calculation uses (1-1/N)[mu2_lower-h(mu2_upper+1/25)] and mu2_upper+1/(250Nqmin). Section 8 proves their whole-tail validity. Mean distance to either excluded endpoint in event 6 is at least x/200000. The independent free target coordinates and J_2<1/(125rho_b) give variance <=n/(125^2rho_b^2). Thus failure is

    <=5760000/n<=5760000/N<=epsA/8.                    (14)

The six failures total at most 6epsA/8<epsA. Therefore

    lambda(A)>=(1-epsA)n log(rhi/rlo)>mL n.             (15)

These estimates precede conditioning on A or D. They impose no ambient distribution on mu_a.

## 5. Whole matching profiles and true source-face variance

The multiplicative count sign is NEGATIVE:

    ellpin=(1+a0/R)/2,
    ellpin'=-(1-epsilon)h/(2R^2)<0.

Nsum is nondecreasing, so Nsum/(2n ellpin) is nondecreasing. At q_- its mean shortfall from one is <=2h+1/n; at q_+, R=r, its mean upper ratio is <=1. Count fluctuations change the ratio by <=eN/4. The NEW gate 2h+1/N<=3eN/4 hence proves over the complete interval

    1-eN<=Nsum/(2n ellpin)<=1+eN.                       (16)

Using the strong 3/rho derivative from (8), n/rho<=3L/2 and interval width 4h/3 gives log-product shift <=6Lh=6/768. With the baseline event,

    P0<=product_i J_i<=1/P0.

The new gate x/(3qmin N)<=eJ gives all J_i in [1-eJ,1+eJ] pointwise, without a union bound over n. Summing EVERY captured signed source face then gives

    U=[product_i J_i/(2ellpin)]sum_i N_i/J_i,
    Ulo n<=U<=Uhi n.                                   (17)

Both signs and all omitted-coordinate terms are retained.

The remaining strong derivative bounds in (8) give

    |sum J_1(q)-sum J_1(qb)|<=xh,
    |sum J_2(q)-sum J_2(qb)|<=xh/5.

The new gates h<=1/100000 and h/5<1/1000000 imply for all matching q and every w in A

    sum J_1<=2211x/50000,
    5269x/1000000<=sum J_2<=5301x/1000000.               (18)

Fix a captured signed SOURCE face j, before source-cap and soft conditioning. Its pin phi is the deterministic Qpin<=Q0<183/1000. Its other coordinates are independent in the inverse-capture law with EXACT means J_1i/J_i and raw second moments J_2i/J_i. Hence, for V=-log G,

    E(V|j)=Qpin+sum_(i!=j)J_1i/J_i<=p+m_*x,
    Var(V|j)=sum_(i!=j)[J_2i/J_i-(J_1i/J_i)^2]<=v_*x,    (19)
    p=183/1000, m_*=2211/[50000(1-eJ)], v_*=531/100000.

For the upper variance, add the omitted nonnegative J_2 term back and use (18) and the new exact gate (5301/1000000)/(1-eJ)<v_*. This controls the TRUE centered variance. It does not add a pin variance or treat J_2 as the whole V^2. No lower variance or normal approximation is needed.

## 6. Finite per-face Cantelli and the whole-tail retained fraction

Let x0=(231/256)8192=7392 and a=chi-m_*. The exact new base certificate proves a>0, u0=a x0-p>0 and

    Ret=9717/10000<u0^2/[v_*x0+u0^2].                  (20)

For x>=x0, u(x)=a x-p>0 and

    d[(a x-p)^2/x]/dx=a^2-p^2/x^2>0.

Thus u(x)^2/[v_*x+u(x)^2] is at least its base value. Apply Cantelli separately in EACH captured signed source face using (19):

    P_inv(V>chi x|j)<=v_*x/[v_*x+u(x)^2]<=1-Ret.        (21)

Only then mix these uniformly bounded probabilities over faces. Their means need not agree; mixture variance is not identified with within-face variance.

The ratio of bsoft to theta/G is min(1,G/theta), which equals one on G>=theta and lies in [0,1] elsewhere. The COMPLETE soft source-shell coefficient therefore lies between Ret theta Ulo n and theta Uhi n. All this is under the fullsource inverse law. The original source cap is paid separately next.

## 7. Complete degrees, distinct exponent rates, and finite payments

The exact matching source log-width is log((r+h)/(r-h)). The target factor arad=r/rhi and 1<=atanh(u)/u<=1/(1-u^2), u=h/r, with the new radial-width gate give

    (1-eRad)2h/rhi<=(r/rhi)log((r+h)/(r-h))
                                   <=(1+eRad)2h/rhi.

Thus raw rows lie between 2Ret(1-eRad)Ulo k and Mk. For the column, fix a source and matched r>=R. If ell_i(r)<1 exactly one pinned TARGET sign is captured; otherwise there are at most two. Since ell_i(r)>1/2, the complete target-face coefficient is at most

    (1/2)sum_i T_i product_(j!=i)ell_j(r)
        <=n product_i ell_i(r)<=nG(z).

Multiplying by arad cancels dr/r exactly to dr/rhi. The ordinary matched target width is <=2h, so bsoft G<=theta gives the COMPLETE column bound C_lambda<=2nh theta/rhi=2k. Inactive columns are zero. Thus M>2 controls both ambient degrees after all restrictions.

For this NEW source band set

    eta=(1-qmax)/(2qmax)=1231/23138,
    C_log=log(qmax/qmin).

The original source-cap complement excludes negative pinned faces and has free fraction 1-(1-q)/(2rho). Its full source mass is at most

    C_log n exp[-(1-1/n)eta x].

Both kernel weights are <=1, so this UNWEIGHTED mass pays its deletion from every raw row. The NEW distinct exact exponent rates are

    cS=(eta-chi)(1-sigma)=13627383/3790929920>0,
    cT=sigma/2-chi(1-sigma)=1447/327680>0.               (22)

The new finite gate eta*x/N<=1/1025 and exp(u)<=1/(1-u) give

    exp(eta*x/n)<=1025/1024.

Let C_log^+ be the new exact upper logarithm enclosure from the verifier. Retaining this logarithmic band mass, the cap loss divided by k is at most

    (1025/1024)C_log^+ rhi 768L exp(-cS L)
       <=(1025/1024)C_log^+ rhi 768L/S64(cS L)<=eS,     (23)
    S64(y)=sum_(j=0)^64 y^j/j!.

This uses the FULL positive prefix, not a single Taylor term or a decimal exponential.

For the one original threshold, the accepted cap identity and Icap<=log t give B_n>=(n/2)log t. For L>=8192, exp(L/2)>=(L/2)^2/2>=L>=x, so log t>=L/2 and B_n>=nL/4>gamma n. Also B_n<=nL and A_n<=exp x=n exp(-sigma L). Set

    tau=sqrt[A_n(B_n-gamma n)].

It is positive and has the NEW finite payment

    tau/k<=rhi 768L^(3/2)exp(-cT L)
       <=rhi 768L^(3/2)/S64(cT L)<=eT.                 (24)

The base check uses sqrt(8192)<91 ONLY to bound the true base expression. Section 8 proves the full tail. After source deletion, final row minus this SAME global tau is at least b k. Final rows and columns are at most Mk.

## 8. All-tail proof of the 53 newly calculated base gates

The verifier uses N=L^64/64! and recalculates every gate at this tuple. It checks (5), all probability/profile bounds, (20), the new mass/deficit (2), both DISTINCT prefix payments (23)--(24), and the scalar price/remainder in Section 10. No searched model diagnostic is substituted for a gate.

After N is inserted, all geometry and target-probability upper envelopes are constants or positive sums of negative powers of L. In particular x/N is proportional to L^-63; epsilon qmax/h is bounded by a multiple of L^-128; ML,VL have powers L^-62,L^-63; count/profile bad probabilities have power L^-64. Source-radius lower envelopes Nqmin/x grow as L^63. Radial and count-calibration errors decrease. The H2 lower mean envelope in (13) increases and its upper decreases; the lower factors are positive at the base. The log-product loss 6Lh=6/768 is constant. Cantelli retention has its explicit monotonicity proof in Section 6. Thus these gates hold on the whole tail.

For the delicate exponential payments, for y>=2,

    yS64'(y)-S64(y)>=-1+y^2/2>=0,

and for y>=4,

    yS64'(y)-(3/2)S64(y)>=-3/2-y/2+y^2/4>=0.

Every omitted polynomial coefficient is positive. Thus L/S64(cL) and L^(3/2)/S64(cL) decrease on these domains. The NEW base gates cS L0>4 and cT L0>4 are checked separately. Hence (23)--(24) hold for all L>=8192. The finite-n correction bound persists because eta*x/N decreases. No dimension-tail sampling is used.

## 9. Same original finite history, genuine deficit, and all positive radii

Use the unchanged largest-eligible-radius algorithm and true disjoint increments. For every algorithmic dilation mathfrak a>=8n^4, normalize a positive load-support point by its last minimum positive-contribution radius and reflect the last center into the positive orthant. Write alpha=mathfrak a/(mathfrak a+1), c_cap=(mathfrak a beta-1)/(mathfrak a+1), and u=last scaled center.

The genuine uncaptured-increment measure mu_a is dominated by the FULL finite corner envelope lambda_a. Its total mass is the actual uncaptured weight; lambda_a has mass alpha^-n J_mathfrak_a(t). On A and P the envelope lambda_a is EXACTLY alpha^-n lambda, since D lies inside the finite cap by

    alpha-c_cap m-(1-beta m)
        =[(beta+1)m-1]/(mathfrak a+1)>0.

An actual target item j meeting A satisfies c_cap r_j<=||u-w||<alpha+rhi<3. Since mathfrak a>2t gives c_cap>1, r_j<3. An actual source increment has r_i>=rho>4; it precedes the SAME target j. The displacement margin

    alpha beta-kappa-2/(mathfrak a+1)
       =zeta-(beta+2)/(mathfrak a+1)>0

holds for all mathfrak a>=8n^4. Included pairs really capture that target center in its own original eligibility test. Integrating only included pieces costs at most their genuine item weights. Both arad and bsoft are <=1, so the actual weighted row is <=arad(w)tau<=tau on mu_a-a.e. target. No ambient law is imposed on mu_a.

Put c_a=alpha^-n>=1. Energy scales by c_a^2, degrees by c_a, and target mass by c_a. The kernel is supported on A times P with A intersect P empty. In the energy-difference identity for nu_a=lambda_a-mu_a, the two terms are bounded using the ambient degrees and the missing masses nu_a(A),nu_a(P). Their sum is at most the TOTAL full finite-cap deficit Delta_a. Consequently

    Delta_a >=[c_a^2 E-tau c_a lambda(A)]/[c_a Mk]
        >=[E-tau lambda(A)]/[Mk]
        >=lambda(A)b/M>=Gamma n>gamma n.               (25)

This holds for EVERY original finite history, uniformly in algorithmic a. The two weights are post-selection energy tools; true increments, original weights and choices are unchanged. No independent histories or selector convergence are presumed.

For arbitrary f in L1 use |f|. Every point of its strict maximal level set has an actual positive witness radius; its open a-shrunk cube covers that point, and finitely many cover each compact subset. Selected dilates remain true witnesses. The normalization in (25) works for every positive terminal radius. Inner regularity gives the full level set without compact-support or scale-truncation assumptions.

## 10. Complete cover/root cost and the new exact finite price

The original captured budget is tau and the terminal self-weight is <=mathfrak a^-n. Insert (25) into the accepted finite cover/load/witness formula, retaining its ENTIRE cover factor:

    C_n<=[1+(beta+1/mathfrak a)^n/tau]
       [tau+mathfrak a^-n+alpha^-n J_mathfrak_a(t)-gamma n].

For fixed n let mathfrak a tend to infinity ONLY in this family of numeric bounds already uniform over all original finite families. At the chosen tau it gives exactly

    C_n<=(1+A_n/tau)(tau+B_n-gamma n)
          =F(A_n,B_n-gamma n).                         (26)

The original scalar estimates give, with T=n log n-n log log n,

    B_n<=T-n log(1-sigma)+1,
    A_n+2sqrt(A_n B_n)
        <=n exp(-sigma L)+2n sqrt L exp(-sigma L/2).

The NEW exact price bound is -log(1-sigma)<257/2500. It is not the coarser sigma/(1-sigma), or the old tuple's price. By (3), after paying c=47/250 the scalar allowance is exactly 1/5000. The new base certificate proves

    1/N(L0)+1/S64(sigma L0)+182/S64(sigma L0/2)<1/5000. (27)

This bounds the TRUE remainder 1/n+exp(-sigma L)+2sqrt L exp(-sigma L/2) at the base. The prefix upper bounds decrease throughout the tail: N grows, and the same polynomial inequality with sigma L0/2>=4 also covers powers zero and one-half. Thus the complete cost in (26) is strictly less than T-(47/250)n for ALL L>=8192. This proves the proposed (1), subject to independent acceptance of the full chain.

## 11. Actual run, frozen local sources, and status

`verify.py` was actually run once for this registered NEW tuple and returned PASS on 53 exact rational gates. It uses Fraction arithmetic in ONE Python process, 64 positive atanh-log terms with exact geometric tails, and FULL S64 exponential prefixes. No producer helper, old passed gate batch, grid, quadrature, shared MATLAB, n=7, lower-bound search or Git operation was used.

The actually read frozen local mathematical sources are:

* R064 relative J/AC/baseline lemma: `c60d48d116863c9c9cac177e19a9fa4158f12fd8dcef0d5239a5ad0a5487d9f1`.
* R066 first-profile lemma: `67bdf01d9aa65c70bf058df08b527b904e954d19f3327ecbb43e5b5c316b5030`.
* R070 second-profile LOCAL: `4d2f4fc2d0f3353ba82e32487dba68804d3021d4c0d477137c11f6c1a004e2b8`.

Their paths and hashes are checked before computation in the verifier. Their abstract active guards/law identities are checked here; notably the R070 Section-7 h_c interval is NOT used to infer the hden=768 finite domain. These local sources do not themselves approve the calibrated target mass/column, same-history bridge or full operator theorem above.

`registration.json`, `registration_amendment.json` and `verification.json` preserve actual candidate and run identities. Only this owned R073 directory was written. The author does not mark the proposed result accepted before independent review.
