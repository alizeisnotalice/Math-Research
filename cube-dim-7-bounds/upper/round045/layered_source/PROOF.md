# R045 — inverse source-column products and a weighted-row second moment

Status: complete proposed proof, submitted for independent review. The accepted R042/R043 sources are unchanged. The proposed improvement uses source weights greater than one and pays their actual eligibility cost explicitly; it is outside R043's uniform-column-cut template. No production scan, n=7 work, MATLAB run, lower-bound search, or Git operation is performed here.

## 1. Proposed theorem and fixed constants

Write L=log n and T=n log n-n log log n. Fix rational constants

    d=2^(-64),   epsilon=2^(-10),   K=2^20,
    B=4K=2^22,   H_R=2^24.

For every sufficiently large n put

    x=(1-Bd^2)L,  t=n/x,  beta=1+1/t,
    h=epsilon/(xd),  zeta=n^(-3),  kappa=beta-zeta,
    I=[d/4,3d/4],  M=nh,  Q=exp(Kxd^2).

The proposed conclusion is

    C_n <= F(A_n(t),B_n(t)-2^(-32)nd),                     (1)
    C_n <= T-2^(-97)n                                    (2)

for all sufficiently large n. Here F(A,B)=(sqrt A+sqrt B)^2, and A_n,B_n are exactly the old accepted scalar family. No explicit finite threshold is claimed. In particular this note does not reuse R042's certified numerical threshold for its different parameters.

The essential intermediate statement is a weighted energy deficit of size 2^(-32)nd. Source inverse normalization replaces the old mean capture scale nh*exp(-Theta(xd^2)) by M=nh. The source-weight ceiling Q multiplies the true eligibility threshold, but the deliberate scalar buffer B>2K makes Q*tau/M negligible.

All eventual guards below are with these fixed constants. In particular x is a fixed positive fraction of L, t tends to infinity, h<d/16, zeta*t<h/2, kappa>1, and xd>=1 eventually.

## 2. Full target/source measures and the one chosen source statistic

Use the accepted corner-cap set

    D={||z||_infinity<t: exists i, z_i<1-beta*m(z)},
    m(z)=max(1,||z||_infinity).

Let lambda^0 have density 2^(-n)m(z)^(-n) on the relevant full annuli, and lambda=lambda^0 restricted to D. A full shell face at radius R has normalized mass dR/(2R), with its other n-1 coordinates independent Uniform[-R,R]. The complete shell mass is n dR/R.

Target annulus A0 has radius r=1+u, u in I. Source annulus P0 has radius rho=tq, q=1-v, v in I. They are disjoint eventually, with r<2 and rho>4. A target/source pair matches when

    |u+v-d|<h.                                           (3)

For a fixed source radius define

    q0=kappa*rho-rho=q-zeta*rho,
    R=q+d-h,       g=R-q0=d-h+zeta*rho.

Every matched actual target radius satisfies

    R<r<R+2h,   1<R<2,   d/2<g<=d,   q0>1/2             (4)

eventually. Define the normalized target-coordinate capture length

    ell_rho(s)=1-(|s|-(rho-g))_+/(2R),     |s|<=rho.

It equals the fraction of w in [-R,R] with |s-w|<kappa*rho. In particular

    1-d/2 <= ell_rho(s)<=1.

The single source statistic and importance weight are

    Phi(z)=-sum_i log ell_rho(z_i),
    omega(z)=exp(Phi(z))=1/product_i ell_rho(z_i).          (5)

For a fixed n this is a bounded measurable function, although its untruncated maximum can be large. The actual kernel later retains only Phi<=Kxd^2, so its weight is at most Q. The source layers are levels of this explicit coordinate statistic, not layers of an unknown actual-mu density.

For now use the full weighted kernel on A0 x P0,

    K0(w,z)=omega(z)*1_match*1_{||z-w||_infinity<kappa*rho}.

Let R0(w) and C0(z) be its full-lambda^0 row and column degrees.

## 3. Pointwise column normalization

For a fixed source z and a matched target radius r, write a_i=kappa*rho-|z_i|>0. The target coordinate capture fractions are

    ell_i(r)=min(1,(1+a_i/r)/2).

They are nonincreasing in r, and r>=R. All fractions are at least 1/2. The full target shell coefficient is at most 2n times their product: each deleted-coordinate product is at most twice the full product, and the normalized face factor is 1/2. The matched target interval has ordinary width at most 2h and logarithmic width at most 2h/R. Hence

    C0(z) <= omega(z)*(4nh/R)*product_i ell_rho(z_i)
          <=4nh=4M.                                     (6)

This holds for every source, without source concentration and before any source cutoff. It is crucial that deleting a source set after inverse normalization costs its original lambda^0 mass times 4M, not that mass times Q.

## 4. Exact coordinate row law and elementary bounds

For a source radius rho and a target coordinate w in [-r,r] define

    J_rho(w)=E_{Z uniform[-rho,rho]}
       [1_{|Z-w|<rho+q0}/ell_rho(Z)].                    (7)

It is even and nonincreasing in |w|. Its mean for a target uniform on the reference interval [-R,R] is exactly one, by Fubini and the defining target capture fraction ell_rho:

    E_{w uniform[-R,R]} J_rho(w)=1.                      (8)

Here is a direct bound which does not rely on numerical integration. Let F_rho=E[1/ell_rho(Z)]. Since ell>=1/2 and the free-coordinate mean loss is g^2/(4Rrho),

    0<=F_rho-1<=g^2/(2Rrho)<=2d^2/rho.

If |w|<=q0 all sources are captured, so J=F_rho. Otherwise the excluded one-sided source interval has length |w|-q0<=g+2h<=3d. Since 1/ell<=2, its normalized integral is at most 3d/rho. Thus, for every actual target coordinate,

    1-3d/rho <= J_rho(w)<=1+2d^2/rho.                    (9)

For a free target coordinate W uniform[-r,r], expansion from the reference radius R, using monotonicity and (8), gives

    E J_rho(W)=1+(1/r)int_R^r [J_rho(w)-1]dw,
    1-6hd/rho <= E J_rho(W)<=1.                         (10)

The boundary region |W|>q0 has probability (r-q0)/r<=3d. Off that region |J-1|<=2d^2/rho, and on it |J-1|<=3d/rho. Consequently

    E(J_rho(W)-1)^2
       <=4d^4/rho^2+27d^3/rho^2
       <=31d^3/rho^2.                                  (11)

The constant 31 is deliberately retained; finite integration samples are not used to lower it.

For two distinct matched source radii rho,rho', the same target W is shared. Do not treat the two J factors as independent. Their deviations have nonpositive means by (10), and Cauchy--Schwarz with (11) gives the uniform crossmoment

    E[J_rho(W)J_rho'(W)]
       <=1+31d^3/(rho*rho').                            (12)

This is the small row-product cost, of order n d^3/rho^2=O(L^2d^3/n), rather than the old column entropy of order Ld^3.

## 5. Pinned source faces and the complete row second moment

Fix a target w and source radius rho. The pinned source-coordinate factor is

    ell_*=ell_rho(+-rho)=1-g/(2R) in [1/2,1].

For each coordinate j, the number of captured pinned source signs is

    N_j(w)=1+1_{|w_j|<q0},

apart from irrelevant boundary conventions; in all cases 1<=N_j<=2. Summing every source face, the exact coefficient per d rho/rho is

    S_rho(w)=(1/2)sum_j [N_j(w)/ell_*]
                            product_{i!=j}J_rho(w_i).   (13)

Every J is between 1/2 and 2 eventually by (9). Therefore the **entire pinned-face sum**, not merely one product or one face, obeys pointwise

    (n/4) product_i J_rho(w_i)
       <=S_rho(w)<=4n product_i J_rho(w_i).              (14)

This handles correlations between all source-face contributions.

Now condition on any maximal target face at radius r. One target coordinate is pinned at +-r; its J is at least 1/2. The other n-1 target coordinates are independent uniforms on [-r,r]. By (10), n/rho=x/q<=2x and h=epsilon/(xd),

    E_target product_i J_rho(w_i)
       >=(1/2)(1-6hd/rho)^(n-1)
       >=(1/2)exp(-24epsilon)>=1/4                      (15)

eventually. The matching source logarithmic length is at least h and at most 4h. Combining (14)-(15), every conditional target face/radius satisfies

    E_target R0(w)>=nh/16=M/16.                         (16)

For the second moment integrate two source radii before taking the target expectation. The pinned target coordinate contributes J_rho J_rho'<=4. Each free coordinate has crossmoment (12), so

    E_target[product_i J_rho(w_i)J_rho'(w_i)]
       <=4 exp(31n d^3/rho_min^2)<=12                    (17)

eventually, since rho_min>=t/2 and x=O(L). Applying (14) to both whole face sums and their logarithmic intervals gives

    E_target R0(w)^2
       <=(4n)^2*(4h)^2*12=3072 M^2<=4096M^2.           (18)

Neither (17) nor (18) assumes independent source radii, independent face choices, or independent pinned coordinates. Independence is used only among the n-1 free coordinates on one fixed target face, as supplied by the actual full-shell measure.

The full target annulus mass a0 satisfies nd/4<=a0<=nd. Thus full weighted energy and the integrated row second moment obey

    E0=int_A0 R0 dlambda^0 >=ndM/64,
    int_A0 R0^2 dlambda^0 <=4096ndM^2.                  (19)

## 6. High inverse weights: delete original source mass, not tilted mass

At a fixed maximal source face one coordinate is pinned at +-rho, and its contribution to Phi is at most d. For every free coordinate, -log ell<=2(1-ell), so

    E[-log ell]<=g^2/(2Rrho)<=d^2/(2rho).

Therefore, uniformly in source radius and pinned face,

    E Phi<=d+xd^2<=2xd^2                               (20)

eventually. Markov on this nonnegative **logarithmic** statistic gives

    lambda^0_source{Phi>Kxd^2}<=2nd/K,                  (21)

since the full source annulus has mass at most nd. By the pointwise normalized column bound (6), deleting these sources removes weighted energy at most

    4M*(2nd/K)=8ndM/K=2^(-17)ndM.                     (22)

No Q multiplier occurs in (22). Q is used only in the actual eligibility row later. This is distinct from retaining a fixed fraction of old energy with a cap normalized by the old column mean.

## 7. Target row cutoff and both original corner-cap restrictions

Delete targets with R0>H_R M. Equation (19) gives removed weighted energy at most

    (1/(H_R M))int R0^2 dlambda^0
       <=4096ndM/H_R=2^(-12)ndM.                       (23)

We must also impose the original D cap on targets and sources. For target radii r<1+3d/4, its cap-complement face fractions are at most (3/4)^(n-1) eventually. Its lambda^0 mass is therefore at most nd*(3/4)^(n-1). Cauchy--Schwarz with (19) bounds the removed target energy by

    64ndM*(3/4)^((n-1)/2)=o(ndM).                       (24)

For a source radius rho=tq, v=1-q>=d/4, the original cap complement has each free-coordinate fraction 1-v/(2rho); the negative pinned face is excluded. Hence its total source mass is bounded by

    nd*exp(-xd/16)=o(nd).                               (25)

Use (6), not a source-weight supremum, to bound its removed weighted energy by 4M times (25). Both (24) and this source-cap energy are at most 2^(-12)ndM eventually.

Define the final sets and kernel by

    A=A0 intersect D intersect {R0<=H_R M},
    P=P0 intersect D intersect {Phi<=Kxd^2},
    K(w,z)=omega(z)*1_A(w)*1_P(z)*1_match
                        *1_{||z-w||<kappa*rho}.

Union-bounding the four removed energies from (22)-(25),

    E=int K dlambda dlambda
      >=[1/64-3/4096-1/131072]ndM
      >=ndM/128.                                       (26)

Also lambda(A)<=nd, ambient rows are at most H_R M, and ambient columns are at most 4M<=H_R M. All these are pointwise bounds for the actual final weighted kernel. No source or target distribution assumption has been imposed on the unknown actual measure mu.

## 8. Weighted energy lemma and genuine mu<=lambda transfer

The accepted energy argument extends to any nonnegative bounded kernel; it does not require an indicator or K<=1. Indeed for nu=lambda-mu>=0 on disjoint A,P,

    E_lambda,lambda-E_mu,mu
       =int K dnu_A dlambda_P+int K dmu_A dnu_P
       <=D0*(nu(A)+nu(P))<=D0*Delta,

if both ambient degrees are <=D0 and Delta is the total envelope deficit. If the true weighted row is at most Q*tau on mu-a.e. target, then E_mu,mu<=Q*tau*lambda(A). Consequently

    Delta >=[E-Q*tau*lambda(A)]/(H_R M).                 (27)

This is an algebraic measure inequality, not an average-case substitute for actual eligibility.

For the original finite algorithm use its usual large algorithmic dilation a_alg, alpha=a_alg/(a_alg+1), c=(a_alg beta-1)/(a_alg+1), actual uncaptured-increment measure mu_a and corner envelope lambda_a. As in the accepted R033/R035/R042 bridge, mu_a<=lambda_a, and on A,P subset D one has exactly lambda_a=alpha^(-n)lambda. This follows from D subset H_a because

    alpha-c*m-(1-beta*m)=((beta+1)m-1)/(a_alg+1)>0.

Take a_alg>=8n^4. Then kappa+2/(a_alg+1)<alpha beta. A target point w in A belongs to an actual uncaptured item j whose radius satisfies

    c*r_j<=||u-w||<alpha+1+3d/4<3,

so r_j<3 because c>1. Every source z in P belongs to an item i with r_i>=rho>4. The same original largest-radius history therefore gives i<j, and actual center displacements imply capture of y_j by bQ_i exactly as in the accepted bridge. The unweighted actual row is at most tau by that **one target's** true eligibility.

The new source weight satisfies omega<=Q on P. Hence the same actual eligibility gives the required weighted row at most Q*tau, including all radial, row-cut and source-statistic restrictions. No new selection or f-dependent source density is introduced.

Writing c_a=alpha^(-n)>=1, energy scales by c_a^2, ambient degrees by c_a, and target mass by c_a. The finite-a version of (27) is at least

    [c_a^2E-Q*tau*c_a lambda(A)]/(c_a H_R M)
      >=[E-Q*tau lambda(A)]/(H_R M).

If Q*tau<=M/256, equations (26)-(27) therefore give, for every original finite history and every a_alg>=8n^4,

    Delta>=nd/(256H_R)=2^(-32)nd.                        (28)

## 9. Scalar buffer pays the true source-weight cost

Let delta0=2^(-32)nd and choose the actual selection threshold

    tau=sqrt[A_n(t)*(B_n(t)-delta0)].

The old integral bound I_n(t)<=2t/(n-1) implies B_n(t)>=nL/2 eventually, so B_n>delta0. Also A_n<=exp x=n exp(-Bd^2L) and B_n<=nL. Thus

    tau<=n sqrt L exp(-Bd^2L/2),
    Q*tau/M
       <=(sqrt L/h)exp[-(B/2-K)d^2L]
       =(sqrt L/h)exp(-Kd^2L) ->0.                      (29)

This proves the actual weighted-row condition in (28). The factor Q has not been omitted or absorbed into an unjustified mean.

Insert the uniform finite-a deficit into the original finite cover/load formula:

    C_n <=(1+((beta+1/a_alg)^n)/tau)
           [tau+a_alg^(-n)+alpha^(-n)J_a(t)-delta0].

Let a_alg tend to infinity only in these already uniform scalar bounds. This gives (1). The finite compact-level-set argument is the original one: each x has a true witness radius R_x, the small cubes Q(x,R_x/a_alg) cover a compact subset, and every selected a_alg-dilate remains a true witness. Lebesgue inner regularity then handles the full level set. No support restriction on f, scale truncation, selector limit or convergence of histories occurs.

Finally, with Bd^2<1/2,

    B_n(t)<=T- n log(1-Bd^2)+1<=T+2Bd^2 n+1,
    A_n+2sqrt(A_n B_n)=o(n)=o(nd).

For the stated constants 2Bd^2=(2^(-32)d)/512. Therefore (1) is at most T-delta0/2 eventually, proving the proposed conclusion (2): delta0/2=2^(-97)n.

## 10. Proof status, arithmetic checks and boundary of the claim

The analytic proof above supplies: both ambient degree budgets, the crossmoment for distinct source radii, the full pinned-source-face sum, original-mass control of the source cutoff, the finite-a actual-row bound Q*tau, and the original general-L1/all-positive-scales passage. The new order comes from allowing bounded source weights greater than one and choosing the scalar buffer to pay them. Merely lowering the old uniform column threshold would not have this property.

The only local computation attached here checks exact rational constant inequalities and source identities; it is not a production integration batch and does not replace the analytic coordinate proof. The coordinator's separate R046 numerical integration is supportive diagnostics only, with no interval enclosure and no all-parameter theorem inferred from samples. Constants 6 and 31 are retained.

No explicit sufficiently-large-n threshold is certified here. Every use of 'eventually' is for the single fixed rational d and constants in section 1, so finitely many guards have a common threshold. Independent review of the new weighted source kernel and transfer is required before treating (1)-(2) as accepted results.
