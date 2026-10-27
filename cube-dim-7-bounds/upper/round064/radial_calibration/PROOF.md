# R064 candidate — multiplicative matching and bounded target radial calibration

Status: **complete proposed fixed-tuple mechanism, pending independent
whole-chain review**. The proposed conclusions are, for every sufficiently
large integer n,

    C_n<=F(A_n(t),B_n(t)-(61/200)n),
    C_n<n log n-n log log n-n/5,                         (1)

where F(A,B)=(sqrt A+sqrt B)^2 and A_n,B_n are the accepted original
scalar quantities. No explicit finite threshold is certified. A scalar
certificate or finite diagnostic by itself does not establish (1).

The operator throughout is the centered axis-parallel common-radius cube
maximal operator on Lebesgue measure, for arbitrary L1 input and all
positive radii. Both new kernel weights are at most one. The new local
condition is a relative gap g/R<=1/3; g can exceed 1/3, so the old absolute
gap domain is not used.

## 1. Fixed rational tuple and complete degree budget

Fix

    s=3/2,       qmin=2/3,        qmax=909/1000,
    delta=2^-30, rlo=1+delta,    rhi=s qmax-delta=2727/2000-delta,
    sigma=91/1000,              chi=1/20,
    e=1/1024,                  h_c=2^-64.

All eps_J, eps_N, eps_Ret, eps_S, eps_T, eps_D, eps_A equal e.
The baseline log tolerance and the radial-width tolerance also equal e.
Put

    P0=1-e-2^-50,
    Ulo=P0(1-e)/(1+e),
    Uhi=(1+e)/[P0(1-e)],
    M=2(1+e)Uhi,
    b=2(1-eps_Ret)(1-e)Ulo-eps_S-eps_T,
    mL=1239/4000,
    Gamma_L=mL b/M,
    gamma_budget=(1-eps_D)Gamma_L,
    gamma=61/200.

Here P0 is a scalar product lower bound, distinct from the full source
set mathcal P_0 defined below. The exact fixed-tuple certificate proves

    0<P0<1, M>2, b>0,
    (1-eps_A)log(rhi/rlo)>mL,
    0<gamma<gamma_budget<Gamma_L<1.                      (2)

The full scalar coefficients are exactly

    c_raw=gamma-sigma/(1-sigma)=37249/181800,
    c_raw-1/5=889/181800>0.                              (3)

For L=log n define

    x=(1-sigma)L=909L/1000, t=n/x, beta=1+1/t,
    zeta=n^-3, epsilon=zeta t,
    kappa=beta-zeta, h=h_c/L,
    theta=exp(-chi x), k=nh theta/rhi.

In particular 0<x<=L. The algorithmic small-cube dilation is denoted
mathfrak a below and is unrelated to these fixed constants.

## 2. Full bands, final kernel, and all active local guards

Use the original corner envelope

    D={||z||_infinity<t: some z_i<1-beta m(z)},
    m(z)=max(1,||z||_infinity),
    d lambda=2^-n m(z)^-n 1_D dz.

Before cap restriction write lambda0 for the same density on the full
annuli. A signed maximal face at radius R>1 has weight dR/(2R) with
n-1 independent Uniform[-R,R] free coordinates.

Define

    A0={w:r=||w||_infinity in [rlo,rhi]},
    mathcal P_0={z:rho=||z||_infinity, q=rho/t in [qmin,qmax]},
    P=D intersect mathcal P_0.

Matching is multiplicative:

    |r-sq|<h.                                          (4)

For each source in the FULL band mathcal P_0 define

    R=sq-h,
    a0=(kappa-1)rho=(1-epsilon)q,
    g=R-a0=(s-1+epsilon)q-h,
    ell_R(z_i)=min(1,(R+kappa rho-|z_i|)/(2R)),
    G(z)=product_i ell_R(z_i),
    b_soft(z)=min(1,theta/G(z)).                         (5)

These definitions include sources outside D and sources with no matched
target. On this band eventually R>0, a0>0, and G>0, so 0<b_soft<=1.
Extend b_soft=0 outside mathcal P_0 if a globally defined weight is
desired. Define the TARGET radial factor on A0 by

    a_rad(w)=r/rhi,

and extend it by zero off A0. Thus 0<=a_rad<=1 globally.

After the good target set A is defined in Section 4, the FINAL kernel is
exactly

    K(w,z)=1_A(w)1_P(z)1_{|r-sq|<h}
              1_{||z-w||_infinity<kappa rho} a_rad(w)b_soft(z).

Its ambient degrees are R_lambda(w)=integral K(w,z)d lambda(z) and
C_lambda(z)=integral K(w,z)d lambda(w). All auxiliary fields are Borel
measurable. The full source probabilities below are taken on mathcal P_0
BEFORE its cap restriction to P or soft weighting.

Eventually h<min(delta/2,1/1000), 0<=epsilon<1/2, and

    epsilon qmax<=h/s.                                  (6)

For each target r the whole interval

    q_-=(r-h)/s, q_+=(r+h)/s                             (7)

lies in [qmin,qmax]. Its width is 2h/s. Its right endpoint
q_b=q_+ is the baseline, since R(q_b)=r exactly.

A source is ACTIVE when some r in [rlo,rhi] matches it. For such a
source, R>rlo-2h>1. For every matched pair, R in [r-2h,r] and

    g-R/3=epsilon q-2h/3<=0,
    0<g<=R/3.                                          (8)

Positivity follows from q>=2/3 and h<=1/1000. Notice the positive
epsilon q=zeta rho term is retained. Some inactive full-band sources
have R<1; their final columns are zero. They are not used in the local
R>=1 lemma or its probability estimates.

Eventually kappa>1, beta<2, t tends to infinity, and rho>=tqmin>=4rhi.
Thus on the complete active matching interval the local guards are

    q>=2/3, 1<=R<=r<3/2, rho>=4r,
    0<g/R<=1/3, a0>0, h<=1/1000, epsilon q<=h/s.

Both annuli lie inside ||z||<t, and A intersect P is empty: targets
have r<3/2 while sources have rho>4 eventually. Definitions (4)--(8)
hold independently of an unknown actual increment measure.

## 3. The new relative-coordinate chain, not an absolute-gap substitution

For fixed target coordinate |w|<=r set

    J(w;q)=(2rho)^-1 integral_[-rho,rho]
                       1_{|z-w|<kappa rho}/ell_R(z) dz.

At the source reference radius R, ell_R(z) is exactly its target-
coordinate capture fraction. The stable closed form, with

    F(R,u)=-2Rlog(1-u/(2R))-u,
    m=(|w|-a0)_+, v=(g-m)_+,

is

    J=1+[F(R,g)+F(R,v)-m]/(2rho).                        (9)

Since a0>0 only one source endpoint can be lost. Also m<=g+2h and
rho>=4r, so (9) includes the m>g central-segment deletion when |w|>R.
All F arguments lie in [0,g]; no negative log argument is introduced.

For u<=g<=R/3,

    0<=F_u<=1/5, |F_R|<=1/15, 0<=F<=R/30.

Now R'=3/2 and g'=1/2+epsilon<1. Hence
|[F(R,g)]'|<=3/10. Writing B=F(R,g)+F(R,v)-m, one has a.e.
|B'|<=17/10: when m=0 both F terms are the same; when m>0,
v=(R-|w|)_+ has derivative zero or 3/2 and -m'=1-epsilon<=1.
Furthermore

    |B|/q<=R/(15q)+R/(3q)+2h/q<=3/5+3h<1,
    J>=1-m/(2rho)>1/2.

Since rho'=rho/q, this yields

    |partial_q log J|<3/rho<=256/rho.                    (10)

The finitely many positive-part branches give local Lipschitz continuity
and absolute continuity on the closed matching interval. Capture onset
can have a genuine derivative kink; (10) is an a.e./AC statement. This
bound is unconditional in |w|<=r, so no extra K_w target condition is
needed. Also (9) gives uniformly

    |J-1|<=max(F(R,g),m/2)/rho
          <=(r/6+h)/rho<1/(3rho).                       (11)

The independently checked relative-coordinate lemma cited in Section 10
details the same new derivative chain and its domain.

At q_b, R=r and the BASELINE target law is W uniform[-r,r]. Fubini
gives E J=1. The elementary moment calculation only needs R>=1,
g<=R/2 and rho>=4R; it does not require g<=1/3. In brief, with
T=(|W|-(R-g))_+, the representation J-1=(F(R,g)-A(T)/2)/rho has

    F(R,g)<=g^2/(3R), A(T)<=4T/3, E T^2=g^3/(3R).

Consequently E(J-1)^2<=11g^3/(27Rrho^2)<=3g^3/rho^2.
Together with |J-1|<=g/rho<=1/8, Taylor's log bounds and Jensen give

    -6g^3/rho^2<=E log J<=0,
    Var(log J)<=12g^3/rho^2.                            (12)

The pinned target coordinate +-R has a deterministic log J in
[-g/rho,0]. The fullsource law used in Section 6 is a different law.

## 4. Target typicality and the new decreasing-pin count calibration

Fix a full target face and radius r. At q_b the other n-1 coordinates
are independent Uniform[-r,r]. Since g<=r/3<1/2 and rho_b>=tqmin,
(12) gives the uniform bounds

    |E sum_i log J_i(q_b)|
      <=6(1/2)^3 L^2/(qmin^2 n)+(1/2)L/(qmin n) ->0,
    Var(sum_i log J_i(q_b))
      <=12(1/2)^3 L^2/(qmin^2 n) ->0.

Chebyshev therefore gives a uniform o(1) failure probability for the
baseline event |sum_i log J_i(q_b)|<=e.

For a fixed target w the number of captured pinned source signs at
coordinate i is

    N_i(q)=1+1_{|w_i|<a0(q)}, Nsum(q)=sum_i N_i(q),
    ell_pin(q)=(1+a0/R)/2.

Because h<=1/1000 and r>=1, a0(q)<r over matching; the pinned target
coordinate contributes one. Under the full target-face law,

    E Nsum(q_+/-)=n+(n-1)a_+/-/r,
    Var Nsum(q_+/-)<=n/4,
    a_+/-=(1-epsilon)q_+/-.

Impose both endpoint count events

    |Nsum(q_+/-)-E Nsum(q_+/-)|<=eps_N n/4.               (13)

Each has failure probability at most 4/(eps_N^2 n), uniformly in r.

The new derivative has the important NEGATIVE sign

    partial_q ell_pin=-(1-epsilon)h/(2R^2)<0.             (14)

Nsum is nondecreasing, so Nsum/(2n ell_pin) is nondecreasing. Its
minimum occurs at q_- and maximum at q_+. At q_-, R_-=r-2h>=1 and

    1-[1+(1-1/n)a_-/r]/[1+a_-/R_-]
      <=2h+1/n.

Indeed the numerator discrepancy is
2h a_-/[r(r-2h)]+a_-/(nr), divided by a denominator at least one.
The fluctuation in (13) changes the ratio by at most eps_N/4. At
q_+, R=r, the mean ratio is at most one, and its upper fluctuation
is again at most eps_N/4. Thus, when 2h+1/n<=3eps_N/4,

    1-eps_N<=Nsum(q)/(2n ell_pin(q))<=1+eps_N             (15)

on the WHOLE matching interval. No additive-gap endpoint direction
has been transferred into this proof.

The fourth target condition is the original cap D. Since rhi<3/2,
the cap-complement fraction on a positive target face has each free
coordinate fraction 1+1/(2t)-1/(2r)<=3/4 eventually; the negative
pinned face is excluded. Its bad fraction is at most (3/4)^(n-1).

Let A be A0 intersected with D, the baseline event, and both count
events (13). The FOUR full-target-face bad probabilities are uniformly
o(1), hence their union is at most eps_A eventually. Therefore

    lambda(A)>=(1-eps_A)n log(rhi/rlo)>mL n.              (16)

These probabilities describe the ambient full-shell measure, not any
unknown actual history measure.

## 5. Complete inverse row and the logarithmic radial width

By (10), n/rho<=3L/2, and the matching width 2h/s=4h/3,

    |sum_i log J_i(q)-sum_i log J_i(q_b)|
       <=512Lh=2^-55<2^-50.

The baseline event and elementary exponential bounds give

    P0<=product_i J_i(q)<=1/P0.                         (17)

Also (11) gives 1-eps_J<=J_i(q)<=1+eps_J eventually, uniformly in
all coordinates and matching parameters, without a union bound over n.

Summing EVERY captured signed source face, with inverse weight 1/G,
the complete coefficient per d rho/rho is exactly

    U(w,rho)=[product_i J_i/(2ell_pin)] sum_i N_i/J_i.

Using (15) and (17) proves

    Ulo n<=U(w,rho)<=Uhi n.                              (18)

This retains the full pinned-source factor and all sign/omitted-
coordinate terms.

The exact matching SOURCE log-width is

    integral_q_-^q_+ dq/q=log[(r+h)/(r-h)].

Put u=h/r. The positive atanh series gives

    1<=log[(1+u)/(1-u)]/(2u)<=1/(1-u^2).

Since h/rlo tends to zero, the last bound is at most 1+e eventually.
In particular the width, after multiplying the TARGET factor r/rhi,
has the two uniform bounds

    (1-e)2h/rhi <= (r/rhi)log[(r+h)/(r-h)]
                                  <=(1+e)2h/rhi.        (19)

The lower bound could be sharpened to e=0, but the displayed error is
paid in b. The factor s has canceled in the log-width ratio exactly;
its ordinary q-width remains 2h/s.

## 6. New relative inverse moments and per-face soft retention

Let phi=-log ell_R and V=-log G. For a FULLSOURCE uniform coordinate
Z on [-rho,rho], the two tails have total density du/rho. Substitution
v=-log ell_R gives exactly

    E_full(phi/ell_R)=R Q_pin^2/rho,
    E_full(phi^2/ell_R)=2R Q_pin^3/(3rho),
    Q_pin=-log ell_pin<=log(6/5).                       (20)

The last inequality follows from g/R<=1/3. A five-term odd alternating
Taylor upper bound at 1/5 gives

    log(6/5)<Q5=34187/187500<1/5.

The exact Fraction certificate proves

    s Q5^2/(1-eps_J)<999/20000,
    (2s/3)Q5^3/(1-eps_J)<1/150.                         (21)

On a fixed captured signed source face, the pinned phi is the
deterministic Q_pin, and the other coordinates are independent with
inverse-capture densities 1_capture/(ell_R J_i). Because J_i>=1-eps_J,
dropping their capture indicators in the nonnegative integrals (20)
gives, using R/q=s-h/q<=s,

    E_inv(V | face)<=1/5+(999/20000)x,
    Var_inv(V | face)<=x/150.                           (22)

These are WITHIN-FACE bounds. At x>=8000, the pinned mean is at most
x/40000. Since chi-999/20000=1/20000, the threshold chi x exceeds
the upper mean by at least x/40000. Thus on each captured signed face

    P_inv(V>chi x | face)<=40000^2/(150x) ->0.            (23)

Eventually this is at most eps_Ret. Averaging these uniformly bounded
probabilities over the inverse-row face mixture preserves the same
bound. We do not equate mixture variance to a sum of within-face
variances, and do not assume different face means are identical.

Before cap restriction and before soft weighting, G>=theta on the
retained event. The ratio of b_soft to theta/G is min(1,G/theta).
It equals one on that event and lies in [0,1] everywhere. Hence the
complete SOFT shell coefficient lies between

    (1-eps_Ret)theta Ulo n and theta Uhi n.

Multiplying by r/rhi and integrating (19) proves full-source raw rows
bounded by

    2(1-eps_Ret)(1-e)Ulo k <= raw radial soft row
                              <=2(1+e)Uhi k=Mk.         (24)

All source independence in this section belongs to the fullsource,
fixed-face, capture-and-inverse law. The original source cap is paid
separately next; no cap-conditioned independence is assumed.

## 7. Source cap loss and the exact complete column cancellation

The full source band mass is n log(qmax/qmin)<=Csrc n, where

    Csrc=(qmax-qmin)/qmin=727/2000.

For rho=tq put v=1-q. In the original corner-cap complement, the
negative pinned source face is excluded, and each free coordinate has
fraction 1-v/(2rho). Therefore its mass is at most

    Csrc n exp[-(1-1/n)eta x],
    eta=(1-qmax)/(2qmax)=91/1818.

The exact exponent gap is

    eta-chi=1/18180>0.

Since both a_rad and b_soft are at most one, deletion of this source
cap complement costs at most that UNWEIGHTED mass in any row. Its
ratio to k is bounded by

    (Csrc rhi/h) exp[-((1-1/n)eta-chi)x] ->0.             (25)

For example eventually the exponent coefficient is at least 1/36360;
h^-1 grows only linearly in L. Thus the source loss is at most eps_S k,
and the final ambient rows obey

    R_lambda(w)>=[2(1-eps_Ret)(1-e)Ulo-eps_S]k,
    R_lambda(w)<=Mk,                  w in A.            (26)

For the column, fix a source z and a matched target radius r>=R. Let
ell_i(r) be its target-coordinate capture fraction and T_i(r) its
number of captured pinned TARGET signs. Because kappa rho-|z_i|>0,
ell_i(r)>1/2. If ell_i(r)<1 exactly one target endpoint is lost and
T_i=1; if ell_i=1 then T_i<=2. Thus T_i/(2ell_i)<=1 in every case,
including boundary cases in the upper-bound direction. The COMPLETE
target-face coefficient is at most

    (1/2)sum_i T_i product_{j!=i}ell_j(r)
       <=n product_i ell_i(r)<=n G(z).                  (27)

The last inequality is monotonicity in r>=R. Now the target radial
weight performs the exact cancellation

    (dr/r) a_rad(r)=dr/rhi.

The ordinary matched target interval has width at most 2h. It follows
from (27), with ALL target faces included, that

    C_raw,radial(z)<=2nh G(z)/rhi,
    b_soft(z) C_raw,radial(z)<=2nh theta/rhi=2k.          (28)

An inactive source has zero column. The final good-target and original
cap restrictions only decrease this full column. Thus ALL final ambient
rows and columns are at most Mk, because M>2. The cancellation alone is
not used as a substitute for (27) or for the active reference-radius
capture guards.

## 8. Genuine same-history finite-a energy deficit

Choose the actual original selection threshold

    tau=sqrt[A_n(t)(B_n(t)-gamma n)].

The original integral bound gives B_n>=nL/2 eventually and gamma<1,
so tau is positive. Also A_n<=exp x=n exp(-sigma L) and B_n<=nL, hence

    tau/k <=(rhi sqrt L/h)
       exp[-(sigma/2-chi(1-sigma))L] ->0,
    sigma/2-chi(1-sigma)=1/20000>0.                     (29)

Thus tau<=eps_T k eventually. By (16), (26), and the two ambient degree
bounds, the bipartite energy lemma supplies a deficit at least
Gamma_L n provided the true weighted target row is at most tau. We
establish this requirement for every original finite history below.

Take the algorithmic dilation mathfrak a>=8n^4. At any positive original
load-support point, normalize its last minimum positive-contribution
cube radius to one and reflect the last center into the positive orthant.
Write

    alpha=mathfrak a/(mathfrak a+1),
    c_cap=(mathfrak a beta-1)/(mathfrak a+1),
    u=last center/(mathfrak a+1).

The genuine uncaptured increments give the accepted measure
mu_a=alpha^-n 2^-n sum_i r_i^-n 1_{I_i/(mathfrak a+1)} dz,
dominated by the original corner envelope lambda_a. On A and P,
lambda_a=alpha^-n lambda EXACTLY, since D is contained in H_a by

    alpha-c_cap m-(1-beta m)
       =[(beta+1)m-1]/(mathfrak a+1)>0.

If a target point w in A belongs to an actual uncaptured item j,
c_cap r_j<=||u-w||<alpha+rhi<3. Eventually mathfrak a>2t gives
c_cap>1 and hence r_j<3. Every source point z in P belongs to an
actual item i with r_i>=rho>4. Thus r_i>r_j and i<j in the SAME
original largest-radius history.

The finite-dilation displacement guard is positive for every
mathfrak a>=8n^4:

    alpha beta-kappa-2/(mathfrak a+1)
       =zeta-(beta+2)/(mathfrak a+1)>0.

If K(w,z)>0, the two actual scaled centers differ by less than

    kappa rho+(r_i+r_j)/(mathfrak a+1)
       <[kappa+2/(mathfrak a+1)]r_i<alpha beta r_i.

Hence source i really captures that one target center in its true
eligibility test. Integrating only its included increment pieces costs
no more than its original weight. The two factors a_rad(w),b_soft(z)
are both at most one, so the actual WEIGHTED row remains at most tau
on mu_a-a.e. target. No model distribution for mu_a is assumed.

For completeness, the energy identity with nu=lambda_a-mu_a is
E_lambda_a,lambda_a-E_mu_a,mu_a=E_nu,lambda_a+E_mu_a,nu.
Each term is controlled by an ambient degree, and A intersect P is
empty, so nu(A)+nu(P) is at most the TOTAL corner deficit Delta_a.
Let c_a=alpha^-n>=1. Energy scales by c_a^2, degrees by c_a, and
target mass by c_a. Therefore

    Delta_a >=[c_a^2 E-tau c_a lambda(A)]/[c_a Mk]
       >=[E-tau lambda(A)]/[Mk]
       >=lambda(A)b/M
       >=Gamma_L n>gamma n.                            (30)

This holds for every original finite history and ALL
mathfrak a>=8n^4 once n meets the fixed-tuple eventual guards. It uses
actual increments and one target's real eligibility, rather than
separately optimized histories or independent witness assignments.

## 9. Original cover, all positive scales, and the complete scalar cost

The original captured budget is still tau and the terminal self-weight
is at most mathfrak a^-n. Insert (30) into the accepted original
finite cover/load/witness formula, retaining its cover factor:

    C_n <=[1+(beta+1/mathfrak a)^n/tau]
      [tau+mathfrak a^-n+alpha^-n J_mathfrak_a(t)-gamma n]. (31)

For arbitrary f in L1 use |f|. Every point of its strict maximal level
set has a true positive witness radius R_x. Its open small cube
Q(x,R_x/mathfrak a) covers that point, so finitely many cover every
compact subset. All selected mathfrak a-dilates remain true witnesses.
The normalization leading to (30) is valid for every positive terminal
radius, hence for these unrestricted finite positive-radius families.
Lebesgue inner regularity gives the full strict level set. No compact
support, radius cutoff, fixed scale lattice, selector convergence or
finite numerical model replaces general L1/all-positive-scale input.

For each fixed n, let mathfrak a tend to infinity only in the family
of already uniform NUMERIC bounds (31). The result at the chosen tau is

    C_n<=(1+A_n/tau)(tau+B_n-gamma n)
          =F(A_n,B_n-gamma n).                          (32)

Thus the first proposed conclusion (1) includes the complete cover cost.

Since x=(1-sigma)L, the original scalar estimates also give

    B_n<=T-n log(1-sigma)+1
         <=T+[sigma/(1-sigma)]n+1,
    A_n+2sqrt(A_n B_n)
         <=n exp(-sigma L)+2n sqrt L exp(-sigma L/2)=o(n).

Using (3), the full expression (32) is at most

    T-c_raw n+n[1/n+exp(-sigma L)+2sqrt L exp(-sigma L/2)].

The bracket tends to zero and is eventually less than the exact slack
889/181800. Therefore the proposed second conclusion C_n<T-n/5 follows.
The endpoint c=c_raw is not asserted.

## 10. Eventual status, actual checks, and frozen local source identity

Every parameter is fixed independently of n. All guards above are
uniform over the full radial intervals and all faces. Indeed h tends
to zero, epsilon qmax/h=O(n^-2), t tends to infinity, baseline means
and variances are O(L^2/n), both count failure probabilities are O(1/n),
the inverse retention failure is O(1/x), source cap/k is a polynomial
in L times exp(-c x), and tau/k is a polynomial times exp(-L/20000).
Thus the finitely many eventual guards share one threshold. This note
does not certify an explicit such threshold or transfer a threshold
from an earlier tuple.

The new relative-coordinate lemma was actually read in full at
`/Users/zhengzhihao/Desktop/T/cube-r030-roots-audit/cube-dim-7-bounds/audit/round064/relative_coordinate/LEMMA.md`.
Its frozen SHA256 is
`c60d48d116863c9c9cac177e19a9fa4158f12fd8dcef0d5239a5ad0a5487d9f1`.
It independently verifies the new R'=3/2, g'=1/2+epsilon chain, the
negative ell_pin derivative, AC interface, relative baseline moments,
and fixed-source-face inverse law. It does not certify the target radial
column/count or full operator conclusions; those responsibilities remain
in Sections 4--9 above and require whole-chain independent review.

Run `python3 certify.py` in this directory for the ACTUAL new one-tuple
Fraction check. It uses 64 positive atanh terms plus a rational geometric
tail for the target mass log, and the explicit odd-five Taylor pin-log
upper bound. It verifies (2), (3), (21), and both exponent gaps exactly.
It checks the frozen local-lemma hash before computation and refuses a
different dependency identity. This is a scalar/local-rate certificate,
not a geometric or operator proof by numerical sampling.

No old grid, integration batch, production scan, shared MATLAB, n=7
work, lower-bound search, Git operation, or change to earlier rounds was
performed here. Only this owned R064 directory was written. The proposed
full result (1) is not marked accepted by its author.
