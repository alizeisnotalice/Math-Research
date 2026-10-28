# R080 full weighted-cap theorem — independent ACCEPT

Status: **ACCEPT** for the fixed branch

    C_n <= Phi80=F(A_n(t),B_n(t)-(39029/100000)n)
        < n log n-n log log n-(5473/20000)n,
    t=n/[(89/100)log n], log n>=4096.

The operator is the original centered axis-parallel Lebesgue cube
maximal operator, for arbitrary L1 input, every strict positive level,
and every positive witness radius. This acceptance concerns the complete
new weighted source-cap argument and its original-operator transfer.
It does not rely on an unweighted-cap payment or require the obsolete
unweighted rate to be positive.

All source files were read only. This is an independent temporary
staging directory, **not a Git clone**. New files are confined to
`audit/round080/full`. No original source or receipt was modified,
no commit or push was made, and no shared MATLAB was used. The existing
53- and 18-gate programs were not run. Their receipts were read as
finite arithmetic evidence, not treated as proofs of transfer.

## 1. Read set, locked tuple and actual independent arithmetic

Read in full: the proposed R080 PROOF.md, verify.py, REGISTRATION.json,
and its verification.json; the adjacent R079 weighted-cap proof; and
the selected producer case_988_1000.json. The original R069 threshold
proof was additionally read to confirm the exact old comparison tuple
sigma69=83/1000, gamma69=317/1000 and domain L>=2^30. R073 and R078
full proofs and the accepted R073 independent review were already read
in the preceding threshold audit. They supply the accepted local
relative-gap and finite-cover framework, not new numerical thresholds.

The operative R080 tuple is

    qmax=247/250, qmin=2/3, s=3/2,
    sigma=11/100, chi=51/1000, h=1/(1536L),
    rlo=1+2^-19, rhi=741/500-2^-19,
    Ret=999999/1000000, mL=196691/500000,
    gamma=39029/100000, price=5827/50000,
    c=5473/20000, gamma-price-c=1/10000.

The eB/eN/eT=2^-19 and eJ/eRad/eS/eD=2^-20 budgets and epsA=2^-16
are those explicitly restated in the proof. The selected case's old
scope string mentioning R079 and its negative `cap_rate` are stale
metadata; the active weighted-cap exponent, new gamma, target mass,
selected c and complete slack agree with this registered theorem. The
selected receipt hash agrees with REGISTRATION.json, and the compact
receipt registration hash matches its current unchanged registration.

New targeted local Python Fraction arithmetic actually ran, PID 28165,
using an independent N16=L^16/16!, 20-term atanh logarithm enclosure
and P40 positive exponential prefix. It did not repeat all source gates.
The sensitive new cap/mass/deficit/scalar and old-cost separator groups
passed. Exact fractions and all source hashes are in sensitive_result.json.
Display-only values are

    nu=0.6666634877522787,
    K0>14.75878194593,
    weighted cap/k<7.817622887*10^-7<2^-20,
    weighted cap budget margin>1.719120277*10^-7,
    paid Gamma-gamma>1.890997165*10^-5,
    target mass margin>1.329871535*10^-6,
    complete scalar remainder<1.003*10^-44<10^-4,
    Phi69 separating coefficient>0.03964999627,
    Phi78 separating coefficient=0.07131796875.

The N16 check is only additional sensitive arithmetic. The full
all-tail proof below audits the source's Nstar=N32 envelopes, whose
fresh compact receipt is preserved without replay.

## 2. Complete active geometry and pre-cap laws

The change in qmax raises rhi but keeps rhi<3/2. At the new base the
proof's independent compact receipt checks all the active domain guards,
including h<=2^-20, 2h<delta, epsilon*qmax<=h/s, rho>4rhi, x/n<1,
and x*n^2>1. With q_±=(r±h)/s, q_b=q_+, the COMPLETE matching
interval lies in [qmin,qmax] because s*qmin=1 and both endpoint trims
exceed h. At every matching source,

    R=sq-h in [r-2h,r], R>1,
    g=R-(1-epsilon)q>0,
    g-R/3=epsilon*q-2h/3<=0.

The R073 relative-gap local lemmas therefore apply with their actual
strong AC derivative constants, including at positive-part kinks.
Inactive sources may have R<1 but have zero final columns, so they
are not improperly covered by an R>=1 local assertion. On the whole
source band R>0, a0>0, G>0, q<1 and rho<t. A and P are disjoint
because targets have norm<3/2 and sources norm>4.

At baseline R=r, the full-target Fubini identities still hold exactly.
The actual six conditions are baseline log, two endpoint counts,
target D, first profile, and second profile. Their elementary variance
and pin bounds are checked for the new tuple in the compact receipt.
They are all statements on full target faces before conditioning;
independence is needed only among free coordinates on one face.
Six failures cost <=6epsA/8, yielding

    lambda(A)>=(1-epsA)n log(rhi/rlo)>mL*n.

The finite-h second-profile mean uses the exact relation
a-5/6=h/(3R)-epsilon*q/(2R), not only g/R<=1/3.
Endpoint counts and ellpin'<0 give the count ratio on every matching q,
not merely at q_b. AC integration then gives the full product and
profile bridges, with loss 6Lh=6/1536 and errors xh,xh/5. Individual
J_i are separately bounded by 1±eJ. Summing every captured sign and
omitted source coordinate gives Ulo*n<=U<=Uhi*n.

For each captured source face, inverse capture and inverse weight
factor over all its FREE coordinates. The actual normalized density is

    1_capture(z_i)/(2rho*J_i*ell_R(z_i)), i!=j.

The fixed source pin is deterministic, regardless of whether it
coincides with the target pin. Thus the free means J1_i/J_i and raw
second moments J2_i/J_i give the TRUE centered variance bound vx,
not a whole-face second moment or mixture variance. Since ell>=5/6,
the free losses are in [0,log(6/5)] and their centered values have
absolute value <=183/1000. R078's bounded-moment mgf argument applies
uniformly to each such product law. Its base Bernstein exponent is
14.430577..., the S32 base prefix exceeds 10^6, and its displayed
positive derivative extends Ret=0.999999 over all larger x. Only then
are faces mixed. The soft ratio is one on G>=theta and nonnegative
elsewhere, giving the claimed full-source raw row lower coefficient.

## 3. New weighted-cap argument: strip, arbitrary pin, product and dq/q

Write I={i:|w_i|<a0} and mcount=|I|. On the WHOLE matching interval,

    n+mcount>=(1-eN)n(1+a0/R)>=(5/3)(1-eN)n,
    mcount>=nu*n, nu=2/3-(5/3)eN.

Because q<1 and rho<t, the complement of the original corner cap on
the source shell is exactly

    D^c={all z_i>=-rho+d}, d=1-q>0.

Every negative pinned source face has zero D^c mass. Fix ANY captured
positive pinned face j. The full lower strip [-rho,-rho+d) is captured
in each counted FREE coordinate i in I minus {j}:

* The lower capture endpoint is w_i-rho-a0<-rho since w_i<a0.
* The upper endpoint is w_i+rho+a0>-rho+d because
  2rho+w_i+a0>2rho-r>7r>1>d.

These inequalities apply to every source point in the strip, rather
than just one representative or its average. Strict boundary points
have zero mass for the actual absolutely continuous source-coordinate
law. Their use is also harmless for genuine increments dominated by
the ambient envelope.

Its inverse-law strip probability is therefore at least

    d/[2rho(1+eJ)],

using ell<=1 and J_i<=1+eJ. At most ONE of the counted coordinates is
lost to the arbitrary positive source pin. Thus there are at least
mcount-1>=nu*n-1 useful free coordinates, even when j is a counted
coordinate or is unrelated to the target pin. Pre-cap product-law
independence and 1-y<=exp(-y) give

    P_inv(D^c | j)
      <=exp[-(mcount-1)d/(2rho(1+eJ))]
      <=exp[-(nu-1/n)*eta*x/(1+eJ)],
    eta=(1-qmax)/(2qmax)=3/494.

The required positivity nu-1/n>0 is verified using Nstar; without it
one would only have the trivial bound. Other free-coordinate cap
conditions were discarded in an upper-bound direction. Negative
faces have probability zero. The common exponential bound can be
multiplied by each nonnegative inverse face coefficient and summed;
the coefficients need not agree. No independence is asserted after
conditioning on D or after clipping, and none is assigned to mu_a.

The crucial pointwise inequality bsoft<=theta/G holds even on the
clipped portion G<theta. It converts each captured cap-complement
integral to theta times its true pre-cap inverse integral. Summing
ALL source faces gives shell loss <=theta*Uhi*n times the common
cap probability. The shell normalization is the same U appearing
in the raw source row: each signed face has dq/(2q) and U already
includes the one-half in the signed-face sum. There is no missing
factor two or second face normalization.

For each fixed target, integrate over its ACTUAL complete matching
interval, giving log(q_+/q_-)=log((r+h)/(r-h)), rather than a full-band
mass. The radial factor yields

    (r/rhi)log((r+h)/(r-h))<=2(1+eRad)h/rhi.

After division by k=n*h*theta/rhi, theta cancels exactly, and the
complete deletion is at most

    M*exp[-(nu-1/n)*eta*x/(1+eJ)]<=M/S32(K(L)).

This is a bound for the actual weighted capture integral removed by
the cap, not for its unweighted ambient mass. The old cS is negative
(-987633/24700000), which is immaterial to this argument. No old
unweighted Clog/h prefactor, no subtraction by chi in the exponent,
and no unweighted-cap probability is used anywhere in acceptance.

## 4. Complete degrees, tau and every larger L

The raw soft row is bounded below by 2Ret(1-eRad)Ulo*k. The exact
radial-width upper bound gives row <=M*k. For the column, r>=R on
matching target shells and every target pin sign is counted;
T_i/(2ell_i(r))<=1 and product ell_i(r)<=G. The arad factor cancels
dr/r, the matched ordinary width is <=2h, and bsoft*G<=theta. Thus
the entire column is <=2k<M*k, including all target signs; inactive
columns vanish and restrictions only decrease these upper degrees.

The old scalar B lower bound remains B>=n log t/2>=nL/4>gamma*n,
so the SAME numeric threshold tau=sqrt[A(B-gamma*n)] is positive.
The A and B upper bounds imply

    tau/k<=rhi*1536*L^(3/2)*exp(-cT L), cT=961/100000.

The compact source run checks the new S32 bound at L0 using
sqrt(4096)=64, obtaining <eT. This is independent of the weighted
cap argument. Paying the new weighted cap loss <=eS*k and this same
tau leaves row-tau>=b*k with the stated positive b and Gamma=mL*b/M.

All adverse Nstar=N32 geometry/probability envelopes are constants
or negative powers of L: x/Nstar scales as L^-31, epsilon*q/h as
L^-64, baseline log moments as L^-30,L^-31, and bad-event bounds as
L^-32. Source lower radii increase; count, pin, profile-transport and
radial errors decrease. The second-profile lower mean is a product
of positive increasing factors, and its upper mean decreases. This
proves all these base gates on the entire tail.

Bernstein has its own explicit positive derivative. For the new cap,

    K(L)=eta*(89/100)/(1+eJ)*[nu*L-32!*L^-31],

has strictly positive derivative, and K0>4. A positive-coefficient
S32 therefore increases, so the weighted cap bound decreases. There
is no residual polynomial factor. For tau, the polynomial identity

    yS32'(y)-(3/2)S32(y)>=-3/2-y/2+y^2/4>=0, y>=4,

shows that L^(3/2)/S32(cT L) decreases. Its domain holds already at
4096. The same coefficient argument covers the powers zero and
one-half in the complete root remainder. The constant 64 is used
only at the base, never as a uniform upper bound for sqrt L.

## 5. Same finite history, finite a, general L1 and all scales

The new deletion estimate changes no witness, item weight, algorithmic
choice, target set or kernel weight range. With a_alg>=8n^4,
alpha=a_alg/(a_alg+1), ccap=(a_alg*beta-1)/(a_alg+1), the original
finite envelope equals c_a*lambda on the kernel support, where
c_a=alpha^-n>=1. D is inside its finite cap because

    alpha-ccap*m-(1-beta*m)=[(beta+1)m-1]/(a_alg+1)>0.

Since t<n and a_alg>2t, ccap>1. An actual target item meeting A has
radius<3, and an actual source increment in P has radius>=rho>4.
Consequently it precedes that SAME target in the original largest-
eligible-radius history. The displacement guard

    alpha*beta-kappa-2/(a_alg+1)
      =n^-3-(beta+2)/(a_alg+1)>0

holds because beta<2 and a_alg>=8n^4. Every included pair thus
captures the target center in its own original eligibility test.
Integrating only included pieces costs no more than the original
weights; arad,bsoft<=1 give actual row<=tau on mu_a-almost every
target. Auxiliary ambient laws never become actual-increment laws.

Energy scales by c_a^2, degree by c_a and target mass by c_a. For
nu_a=lambda_a-mu_a>=0, the two energy-difference missing-mass terms
cost at most c_a*M*k*[nu_a(A)+nu_a(P)]<=c_a*M*k*Delta_a, since A
and P are disjoint. Actual energy is <=tau*c_a*lambda(A). Therefore

    Delta_a >=[c_a^2 E-tau*c_a*lambda(A)]/(c_a*M*k)
       >=[E-tau*lambda(A)]/(M*k)
       >=lambda(A)*b/M>=Gamma*n>gamma*n.

This holds for every finite history and every such finite dilation.
The c_a factors have the correct direction, and no selector limit is
used. The original cover factor and terminal self-weight remain,
giving the full numeric finite inequality in the author proof. Only
its scalar a_alg limit is taken, yielding exactly F(A,B-gamma*n).

For arbitrary f in L1, |f| has a true positive witness at every point
of each strict level set. The open shrunk cubes give finite subcovers
of compact subsets, and selected dilates remain true witnesses.
Normalization applies to any positive terminal radius. Inner
regularity and level homogeneity pass to every L1 input and every
positive level, without compact support, radius cutoffs, a scale
lattice or convergence of chosen families.

Finally the complete scalar remainder

    1/n+exp(-sigma L)+2sqrt L exp(-sigma L/2)

is paid by the new base prefix check and the whole-tail decrease
above. With -log(1-sigma)<price and gamma-price-c=1/10000, the entire
F value is strictly below T-c*n. The positive A/cross/self terms have
all been included.

## 7. Genuine comparison of complete Phi69 and Phi78

For the original exact B function, the substitution
v=1/rho-1/t gives

    I_n(t)=integral_0^(1-1/t) (1-v/2)^(n-1)/(v+1/t) dv
       <= t*integral_0^infinity exp[-(n-1)v/2] dv
       =2t/(n-1).

Because the core 1-[(1+1/t)/2]^n is nonnegative and n/(n-1)<=2,

    B_n(t_old)>=n log t_old-2n/x_old
       >=T+sigma_old*n-4n/L,

using x_old>=L/2 and -log(1-sigma_old)>=sigma_old. On the accepted
old positive-root domains F(A,B-gamma_old*n)>=B-gamma_old*n. Hence

    Phi69>=T-n[234/1000+4/L], L>=2^30,
    Phi78>=T-n[599/2000-25/256+4/L], L>=8192.

The strict new upper bound Phi80<T-(5473/20000)n lies below both
these old-function lower bounds on their respective domains. The
base separator coefficients above are positive and increase with L.
Thus Phi80<Phi69 on all L>=2^30 and Phi80<Phi78 on all L>=8192.
These lower bounds concern the old upper-bound FORMULAS, not C_n;
there is no reversal or attempt to compare only gamma/display values.
No old branch is assigned a larger validity domain. This acceptance
does not assert global optimality over all accepted branches or kernels.

No mathematical gap was found in the requested new mechanism or its
full-chain transfer. The metadata remnants noted in Section 1 do not
alter the operative registered tuple, new proof or source identities.
