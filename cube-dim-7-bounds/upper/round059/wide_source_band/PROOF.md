# R059 candidate — wide source band, calibrated faces, and a one-tenth scalar margin

Status: **complete proposed mechanism and fixed-tuple certificate, pending
independent full review**. This note does not treat the R058 ideal-model
value as an accepted operator bound. No explicit finite n-domain is
certified. The original operator remains the centered axis-parallel
common-radius cube maximal operator on Lebesgue measure, for arbitrary
L1 input and every positive scale.

The candidate conclusion, if this proof is accepted, is: for all
sufficiently large integer n, with T=n log n-n log log n,

    C_n <= F(A_n(t),B_n(t)-(1637/12500)n),
    C_n < T-n/10,                                         (1)

where F(A,B)=(sqrt A+sqrt B)^2 and A_n,B_n are the original accepted
scalar quantities. The proof uses the entire cover/load cost, not only
the corner budget. The source weights below are bounded by one.

## 1. One fixed rational tuple and its exact certificate

Fix

    d=391/1905,       a_band=283/2032,       b_band=63/64,
    chi=147/10000,    sigma=3/100,
    e=2^-20,         delta=e,              h_c=2^-48,
    eps_J=eps_N=eps_S=eps_T=eps_D=e,
    eps_Ret=1/4096,  eps_A=1/65536,
    baseline log tolerance=e.

Put

    vmin=a_band*d,    vmax=b_band*d,
    qmin=1-vmax,      qmax=1-vmin,
    Rbase=qmin+d,
    rlo=Rbase+delta,  rhi=qmax+d-delta.

For L=log n set

    x=(1-sigma)L,    t=n/x,       beta=1+1/t,
    zeta=n^-3,       kappa=beta-zeta,
    h=h_c/(Ld),      theta=exp(-chi*x),      k=nh theta.

In particular 0<x<=L. The earlier small-d restriction from R051 is NOT
invoked: d<1/4<1/3, and the local coordinate lemma used below has the
independently checked d<=1/3 domain.

Define the following exact rational budgets:

    P0=1-e-1280h_c,
    Ulo=P0(1-eps_N)/(1+eps_J),
    Uhi=(1+eps_N)/[P0(1-eps_J)],
    M=2Uhi/qmin,
    b=(1-eps_Ret)2Ulo/qmax-eps_S-eps_T,
    mL=159477/1000000,
    Gamma_L=mL*b/M,
    gamma_budget=(1-eps_D)Gamma_L,
    gamma=1637/12500.

The attached independent Fraction certificate verifies all denominators
are positive, M>2, b>0, and

    (1-eps_A)log(rhi/rlo)>mL,
    0<gamma<gamma_budget<Gamma_L<1.                       (2)

It uses 64 terms of the positive atanh-log series with a geometric tail
bound, not displayed decimals. It also verifies, with
Qbar=26979/250000,

    -log(1-d/(2Rbase)) < Qbar < 1/9,
    (Rbase/qmin)Qbar^2/(1-eps_J) < 293/20000,
    (2/3)(Rbase/qmin)Qbar^3/(1-eps_J) < 1/900,             (3)
    eta:=vmin/(2qmax),
    eta-chi=499871/37603070000 >1/100000,
    sigma/2-chi(1-sigma)=741/1000000>0.                   (4)

Finally the CHOSEN deficit, slightly below the complete budget, gives

    c_raw=gamma-sigma/(1-sigma)=121289/1212500,
    c_raw-1/10=39/1212500>0.                             (5)

The certificate stores every rational input and derived budget. No
parameter is reconstructed from a decimal model output.

## 2. Ambient measure, matching, and local coordinate guards

Use the original corner-cap domain and measure

    D={||z||_infinity<t: some z_i<1-beta m(z)},
    m(z)=max(1,||z||_infinity),
    d lambda=2^-n m(z)^-n 1_D dz.

Before cap restriction write lambda0 for the same density on the
relevant annuli. A signed maximal face at radius R>1 has weight
dR/(2R), with its other n-1 coordinates independent Uniform[-R,R].

Targets have radius r in [rlo,rhi]. Define the FULL source band and
its final cap restriction explicitly by

    mathcal P_0={z:rho=||z||_infinity, rho/t in [qmin,qmax]},
    P=D intersect mathcal P_0.

The set mathcal P_0 is distinct from the scalar product lower bound P0
in Section 1. Sources have radius rho=tq, q in [qmin,qmax]. Matching means

    |r-q-d|<h.                                          (6)

For fixed source q define

    R(q)=q+d-h,
    a0(q)=(kappa-1)rho=(1-zeta t)q,
    g(q)=R-a0=d-h+zeta t q,
    ell_R(z_i)=min(1,(R+kappa rho-|z_i|)/(2R)),
    G(z)=product_i ell_R(z_i),
    b_soft(z)=min(1,theta/G(z)).                         (7)

Define G and b_soft by (7) on ALL of mathcal P_0, including points
outside D and sources with no matched target radius. They may be extended
by b_soft=0 off mathcal P_0. All these definitions are independent of an
unknown actual increment measure. For each fixed n, G is positive and
the inverse densities used below are finite. On mathcal P_0 the actual
kernel weight is 0<b_soft<=1. The full source laws and inverse-row laws
below are used before the restriction from mathcal P_0 to P.

Eventually h<delta/2, h<d/8 and zeta t<h/2. Hence the FULL interval

    [q_-,q_+]=[r-d-h,r-d+h]                              (8)

is inside [qmin,qmax] for every target radius. Its ordinary width is
2h; its source logarithmic width lies between 2h/qmax and 2h/qmin.
There is no endpoint half-width loss.

If a source has any matched target, then

    q>=rlo-d-h=qmin+delta-h,
    R(q)>=Rbase+delta-2h>=Rbase>1.                       (9)

These are the ACTIVE sources. Inactive source columns are zero. Over
the whole source band R>=Rbase-h>1 eventually, so (7) is well-defined
even for inactive sources. We do not incorrectly apply the active
lower bound Rbase to inactive sources with q very near qmin.

For the active full matching interval, R in [r-2h,r], q>=qmin>1/2,
rho>=4r eventually, a0>0, and

    0<g<=d<R/2.

Thus every guard of the accepted local d<=1/3 coordinate lemma holds.
The source and target annuli are disjoint: r<5/4<2 and rho>4 eventually.
Also t tends to infinity, both annuli lie inside ||z||<t, kappa>1,
and beta<2.

For clarity, the local function is

    J(w_i;q)=(2rho)^-1 integral_[-rho,rho]
                 1_{|z_i-w_i|<kappa rho}/ell_R(z_i) dz_i.

The local lemma supplies on the stated guards:

* At q_b=r-d+h, R=r. Under W uniform[-r,r], E J=1,
  -6g^3/rho^2<=E log J<=0, Var(log J)<=12g^3/rho^2.
  The pinned target coordinate has |log J(+-r)|<=d/rho.
* Uniformly in the full matching interval,
  |J-1|<=5d/(8rho), and log J is absolutely continuous in q.
* Its derivative exists a.e. and is bounded by
  (64/rho)[d^2+1_{|w_i|>r-2d}].

Capture-onset kinks are treated by absolute continuity, not by claiming
classical differentiability at every q.

## 3. Good targets: baseline, the missing binomial event, and two face counts

On a fixed full target face at radius r, the pinned coordinate is +-r
and the other n-1 coordinates are independent uniforms on [-r,r]. At
q_b the log-product mean and variance are bounded uniformly by

    |mean sum_i log J_i(q_b)|
      <=6d^3 L^2/(qmin^2 n)+dL/(qmin n) ->0,
    Var(sum_i log J_i(q_b))
      <=12d^3 L^2/(qmin^2 n) ->0.

Thus the baseline event |sum_i log J_i(q_b)|<=e has failure probability
o(1), uniformly in r and target face, by Chebyshev.

Because d<1/4, the condition K_w<=4nd is NOT automatic. Here

    K_w=#{i:|w_i|>r-2d}
        =1+Bin(n-1,2d/r)

under the same fixed target-face law. Its mean is at most 1+2nd and
its variance at most n/4. For nd>=1,

    P(K_w>4nd)<=1/(4nd^2)=o(1).                         (10)

This is the added event needed to use the whole-interval derivative
bound for the wide-band tuple.

For q in (8) define the number of captured pinned source signs

    N_i(q)=1+1_{|w_i|<a0(q)},
    Nsum(q)=sum_i N_i(q),
    ell_*(q)=1-g/(2R)=(1+a0/R)/2.

Here a0(q)<r. The pinned target coordinate has N_i=1, and each free
target coordinate contributes a Bernoulli(a0/r) in addition to one.
For q_-,q_+ let a_-=a0(q_-), a_+=a0(q_+). Then exactly

    E Nsum(q_+/-)=n+(n-1)a_+/-/r,
    Var Nsum(q_+/-)<=n/4.

Impose both endpoint events

    |Nsum(q_+/-)-E Nsum(q_+/-)|<=eps_N n/4.               (11)

Each has failure probability at most 4/(eps_N^2 n), uniformly in r.
No cap-conditioned coordinate law is substituted for this full-face law.

Both Nsum(q) and ell_*(q) are nondecreasing in q. Consequently (11)
implies on the whole interval

    Nsum(q)/(2n ell_*(q))
      >=1-[2h+1/n+eps_N/4],
    Nsum(q)/(2n ell_*(q))
      <=1+[2hd+eps_N/4].                                (12)

For the lower bound, bound its numerator by Nsum(q_-) and denominator
by 1+a_+/r. The discrepancy a_+-a_- is at most 2h and a_-/(nr)<1/n.
For the upper bound use the opposite endpoints and the exact identity

    a_+/r-a_-/(r-2h)
      =2h(1-zeta t)(d-h)/[r(r-2h)]<=2hd.

Eventually 2h+1/n<=3eps_N/4 and 2hd<=3eps_N/4, so (12) is contained
in [1-eps_N,1+eps_N]. This pays the genuine lost pinned faces rather
than replacing N_i by two.

Finally impose the original target cap D. On each positive target face
its complement has free-coordinate fraction
1+1/(2t)-1/(2r)<=3/4 eventually; the negative target face is excluded
from that complement. Thus the complement fraction is at most
(3/4)^(n-1)=o(1), uniformly in r.

Let A be the target annulus intersected with D, the baseline event,
the K_w event, and both endpoint-count events. The FIVE bad event
probabilities just proved are uniformly o(1). With the fixed eps_A,
their union has probability at most eps_A eventually on every full
target face/radius. Therefore

    lambda(A)>=(1-eps_A)n log(rhi/rlo)>mL n.              (13)

No assumption has been made that an actual history measure has this
random full-face law; (13) is an ambient-envelope statement only.

With this A and the P defined in Section 2, the FINAL kernel is exactly

    K(w,z)=1_A(w)1_P(z)1_{|r-q-d|<h}
                  1_{||z-w||_infinity<kappa rho} b_soft(z),
    r=||w||_infinity, rho=||z||_infinity, q=rho/t.

Write R_lambda(w)=integral K(w,z)d lambda(z) and
C_lambda(z)=integral K(w,z)d lambda(w). These are the final ambient
degrees used below. Eventually A intersect P is empty, because r<=rhi<2
on A and rho>=tqmin>4 on P. A source in mathcal P_0 is inactive when
no r in [rlo,rhi] matches it; its final column is identically zero,
including if it lies in P. This fixes both the pre-cap probability space
and the post-cap energy kernel without changing any estimate.

## 4. Whole-interval inverse row calibration

For w in A, summing the a.e. derivative bound and using K_w<=4nd,
n/rho<=2L gives an upper bound 640Ld. Absolute continuity from q_b
over distance at most 2h hence gives

    |sum_i log J_i(q)-sum_i log J_i(q_b)|
       <=1280Ldh=1280h_c.

The baseline tolerance therefore yields the exact rational enclosure

    P0<=product_i J_i(q)<=1/P0                           (14)

throughout matching, using e^-s>=1-s and e^s<=1/(1-s).
Pointwise |J-1|<=5d/(8rho) and rho>=tqmin tending to infinity also
give, eventually and uniformly,

    1-eps_J<=J_i(q)<=1+eps_J.                            (15)

This pointwise statement needs no union bound over n coordinates.

At fixed source radius, sum ALL captured signed source faces with the
inverse weight 1/G. The exact coefficient per d rho/rho is

    U(w,rho)=[product_i J_i/(2ell_*)] sum_i N_i/J_i.

Combining (12), (14), (15), in their correct lower/upper directions,
proves

    Ulo n<=U(w,rho)<=Uhi n.                              (16)

The pinned source coordinate and its inverse factor are included in
ell_* and N_i. This is the complete face sum, not a one-face bound.

## 5. Fixed-face inverse probability and soft retention

Write phi=-log ell_R and V=-log G. For a uniform free source coordinate
Z on [-rho,rho], the two source tails together have density du/rho,
where ell=1-u/(2R), 0<=u<=g. Substitution v=-log ell gives exactly

    E_full(phi/ell)=R vpin^2/rho,
    E_full(phi^2/ell)=2R vpin^3/(3rho),
    vpin=-log ell_*.

On a given captured signed source face, its n-1 free coordinates are
independent under the normalized inverse-capture law. Each coordinate
density relative to full source law is 1_capture/(ell J_i). By (15)
and nonnegativity, removing the capture indicator for an upper bound
gives mean at most R vpin^2/[rho(1-eps_J)] and second moment at most
2R vpin^3/[3rho(1-eps_J)]. The pinned source coordinate contributes
the deterministic vpin.

For every ACTIVE source, (9), g<=d and q>=qmin imply

    vpin<=-log(1-d/(2Rbase))<Qbar<1/9,
    R/q=1+(d-h)/q<=Rbase/qmin.

The exact bounds (3) consequently prove, on EACH captured signed face,

    E_inv V<=1/9+(293/20000)x,
    Var_inv V<=x/900.                                   (17)

There is no claim that the mixture-of-faces variance equals this
independent-coordinate variance. Instead apply Chebyshev separately on
each face. Since chi-293/20000=1/20000, for x>=40000/9 the distance
from the upper mean in (17) to chi x is at least x/40000. Hence

    P_inv(V>chi x | given captured signed face)
       <=40000^2/(900x) ->0.                            (18)

Eventually this is at most eps_Ret. The same bound then holds for any
inverse-row mixture of these faces by averaging their probabilities,
without assuming their means are identical.

Since theta=exp(-chi x), the event retained in (18) is G>=theta. The
ratio of soft weight b_soft to theta/G is min(1,G/theta), equal to one
on that event and between zero and one elsewhere. Thus the FULL soft
source-shell coefficient is bounded by

    (1-eps_Ret)theta Ulo n <= soft coefficient
                                  <=theta Uhi n.        (19)

All probabilities here are taken before original source-cap or soft
weight conditioning. We next pay the cap deletion separately.

## 6. Source cap deletion and the sharp complete column bound

The full source annulus has mass
n log(qmax/qmin)<=Csrc n, where Csrc=(qmax-qmin)/qmin. For this wide
band one must NOT silently reuse the narrower-band bound nd.

At source rho=tq, v=1-q, the original cap complement excludes the
negative pinned face and retains each free coordinate with fraction
1-v/(2rho). Therefore its full mass is at most

    Csrc n exp[-(1-1/n)eta x],    eta=vmin/(2qmax).

Because b_soft<=1, its removal costs at most this unweighted mass in
ANY target row. Dividing by k gives

    (Csrc/h) exp[-((1-1/n)eta-chi)x] ->0.                 (20)

Indeed eta-chi>1/100000 by (4), so eventually the exponent coefficient
is at least 1/200000, while h^-1 is linear in L. Thus for the fixed
eps_S the cap loss is at most eps_S k. Using the FULL matching width
from (8) in (19), the final ambient rows satisfy

    R_lambda(w)>=[(1-eps_Ret)2Ulo/qmax-eps_S]k,
    R_lambda(w)<= (2Uhi/qmin)k=Mk,        w in A.          (21)

Now fix any source z and an actual matched target radius r>=R. Let
ell_i(r) be its target-coordinate capture fraction. Since
kappa rho-|z_i|>=a0>0, every ell_i(r)>1/2. Let T_i(r) be the number
of captured pinned TARGET signs at coordinate i. If ell_i(r)<1,
exactly one target endpoint is lost, hence T_i=1; if ell_i(r)=1,
T_i<=2. In both cases T_i/(2ell_i(r))<=1. Thus the ENTIRE target-face
coefficient is at most

    (1/2)sum_i T_i product_{j!=i}ell_j(r)
        <=n product_j ell_j(r)
        <=n G(z).

The last inequality is monotonicity in target radius. Its matched
logarithmic target width is at most 2h/R. Consequently

    C_raw(z)<=2nh G(z)/R,
    b_soft(z)C_raw(z)<=2k/R<=2k                           (22)

for an active source; an inactive source column is zero. The original
target/source cap restrictions and good target events only shrink
columns. Since M>2, all final ambient row and column degrees are at
most Mk. The factor two is obtained from the exact pinned target
coupling, not by dropping an uncharged source or target face.

## 7. Tau, energy deficit, and every original finite history

Set

    delta_n=gamma n,
    tau=sqrt[A_n(t)(B_n(t)-gamma n)].

The accepted cap integral bound gives B_n>=nL/2 eventually and gamma<1,
so tau is positive. Also A_n<=exp x=n exp(-sigma L), B_n<=nL, hence

    tau/k <=(sqrt L/h)
       exp[-(sigma/2-chi(1-sigma))L] ->0                 (23)

by the exact positive exponent in (4). Thus tau<=eps_T k eventually.

For a measure mu<=lambda whose actual weighted target rows are at
most tau, the bipartite energy identity on disjoint A and P gives

    Delta >= [E_lambda,lambda-tau lambda(A)]/(Mk),

where Delta is the total envelope deficit. This follows by expanding
lambda-mu in the two variables and using the two ambient degree bounds;
the missing masses on disjoint A and P sum to at most Delta. Equations
(13), (21), (23) imply Delta>=Gamma_L n>gamma n.

We next establish the necessary actual row hypothesis for the SAME
original history and finite algorithmic dilation mathfrak a; the band
parameter a_band never denotes that dilation.

Take mathfrak a>=8n^4. At a positive old load-support point, normalize
the last minimum positive-contribution radius to one and reflect its
center into the positive orthant. Let

    alpha=mathfrak a/(mathfrak a+1),
    c_cap=(mathfrak a beta-1)/(mathfrak a+1),
    u=last center/(mathfrak a+1).

The actual uncaptured increments have their original density
mu_a=alpha^-n 2^-n sum_i r_i^-n 1_{I_i/(mathfrak a+1)} dz.
Their disjointness and center/radius coupling give mu_a<=lambda_a,
the accepted corner envelope. On A and P it is EXACTLY
lambda_a=alpha^-n lambda, since D is contained in H_a by

    alpha-c_cap m-(1-beta m)
       =[(beta+1)m-1]/(mathfrak a+1)>0.

If an actual uncaptured target item j has an increment point w in A,
then c_cap r_j<=||u-w||<alpha+rhi<3. For large n,
mathfrak a>2t so c_cap>1, giving r_j<3. Any actual source increment
point z in P has r_i>=||z||=rho>4, hence r_i>r_j and i<j in the same
original largest-radius history.

The finite-dilation displacement guard is

    alpha beta-kappa-2/(mathfrak a+1)
       =zeta-(beta+2)/(mathfrak a+1)>0.

For each source point captured by the new kernel, its true scaled
center and the target's true scaled center differ by less than

    kappa rho+(r_i+r_j)/(mathfrak a+1)
       <[kappa+2/(mathfrak a+1)]r_i<alpha beta r_i.

So every contributing source item really captures that one target
center in its original eligibility test. Integrating only the included
pieces of its actual increments costs no more than its original weight.
Since b_soft<=1, the new weighted row is at most tau, not merely an
unpaid multiple of tau. This holds on mu_a-a.e. target and uses no
distributional assumption on mu_a.

Write c_a=alpha^-n>=1. Ambient energy scales by c_a^2, degrees by c_a,
and target mass by c_a. The finite-a energy deficit is at least

    [c_a^2 E-tau c_a lambda(A)]/[c_a Mk]
      >=[E-tau lambda(A)]/[Mk]
      >=Gamma_L n>gamma n.                              (24)

Thus every original finite history has the stated corner deficit,
uniformly for ALL mathfrak a>=8n^4 once n is large enough for this one
fixed tuple. No selector limit, new selection rule, independent history
splicing, or convergence of actual increment measures is used.

## 8. Original cover, general L1, all positive scales, and full scalar cost

The captured part of the original load still costs tau and its terminal
self-weight at most mathfrak a^-n. Insert (24) into the original finite
cover/witness argument, retaining its cover factor:

    C_n <=[1+(beta+1/mathfrak a)^n/tau]
      [tau+mathfrak a^-n+alpha^-n J_mathfrak_a(t)-gamma n]. (25)

This original framework applies to every finite family of strict true
witnesses. For arbitrary f in L1, use |f|. Every point of its strict
maximal level set has an actual positive witness radius R_x, and the
open small cubes Q(x,R_x/mathfrak a) finitely cover each compact subset.
Every selected mathfrak a-dilate remains a true witness. The proof of
(24) is invariant under the common normalization of its terminal radius,
so it applies to these arbitrary finite positive radii. Lebesgue inner
regularity then passes to the whole level set. No compact support,
lower or upper scale cutoff, radius lattice, or finite numerical model
is substituted for general L1/all-positive-scale input.

For each fixed n, let mathfrak a tend to infinity ONLY in the numeric
bounds (25) already uniform over every original finite family. Its limit
is

    (1+A_n/tau)(tau+B_n-gamma n)
       =F(A_n,B_n-gamma n)                              (26)

at the specified tau. This proves the first proposed inequality (1),
including the complete cover cost.

Since x=(1-sigma)L, the original scalar estimates give

    B_n<=T-n log(1-sigma)+1
        <=T+[sigma/(1-sigma)]n+1,
    A_n+2sqrt(A_n B_n)
        <=n exp(-sigma L)+2n sqrt L exp(-sigma L/2)=o(n).

Consequently (26) is bounded above by

    T-c_raw n+n[1/n+exp(-sigma L)
                          +2sqrt L exp(-sigma L/2)].

The bracket tends to zero. The exact slack (5) permits it to be less
than 39/1212500 eventually, yielding the second proposed inequality
C_n<T-n/10. A claim of c=c_raw itself is unnecessary and is not made.

## 9. Eventual threshold, provenance, and unresolved review status

All constants in this note are fixed rational numbers. Each guard has
been proved uniformly in target radius, source radius and face. In
particular h tends to zero, zeta t/h=O(n^-2), t tends to infinity,
the baseline means/variances are O(L^2/n), all two-count/K_w failure
probabilities are O(1/n), inverse retention failure is O(1/x), source
cap loss is a polynomial in L times exp(-c x), and tau/k is a polynomial
times exp[-741L/1000000]. Thus the finitely many eventual guards share
one threshold. This note neither computes that threshold nor imports
the explicit R056 threshold for its different tuple.

The current local d<=1/3 J/AC, baseline 3/6/12 moment, and inverse-face
dependency is the final local lemma at
`cube-r030-roots-audit/cube-dim-7-bounds/audit/round056/one_coordinate/LEMMA.md`,
SHA256 `97430af9980c2de5d105bc4ea86667573c5c28e0ff40a92f1c613cd95a243eec`.
That actual final text was read again in full; its local interfaces and
guards were matched to Sections 2--5. `DEPENDENCY_ADAPTATION.md` records
the current binding and the historical receipt identity explicitly.

The unchanged historical `certificate.json` records the source hash
`6d1742b7be914df8e7a26fd0a3f8942ae5275ab0923ce6598d87ad7fd3634800`
observed in its original run. No verified original-text copy of that old
identity was retained, so this revision does not claim that the old source
can be reproduced from the current external path. The receipt is not
rewritten or relabeled as a run against the final lemma. The adaptation
uses the final frozen lemma and does not rerun the 64-term arithmetic.

All new band geometry, target probabilities, face-count calibration,
soft retention, sharp column coefficient and exponent guards were
derived explicitly above. The original finite cover/envelope reference
and the R058 candidate identity remain hashed in the historical receipt.

Actual new computation: `python3 certify.py` in this directory. It uses
local Python standard-library Fraction arithmetic, one fixed tuple,
and rigorous positive atanh-log upper/lower sums. It performs no grid,
quadrature, production scan, MATLAB call or parameter optimization.
All computation is a rational certificate for this proposed proof;
the R058 model alone remains insufficient for an operator conclusion.

All new files are confined to `upper/round059/wide_source_band/` in the
dedicated clone. Old rounds and other agents' files were not changed.
The two claims (1) are submitted for independent review and are not
marked accepted by this author.
