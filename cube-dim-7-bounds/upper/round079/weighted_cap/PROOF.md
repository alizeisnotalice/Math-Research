# R079 LOCAL: weighted source-cap deletion using the existing count gate

Status: proposed local lemma for independent review. It bounds only the
source-cap deletion in the fixed R073/R078 kernel. The finite arithmetic
illustration at L=4096 is conditional on the other active geometry and
profile/count guards holding there; it is not a new complete C_n bound.
The R079 numerical batch using the old unweighted cap model is separate.
No paper, accepted ledger, original file or other agent's file was edited.

Frozen source: R073 `cantelli_frontier/PROOF.md`, SHA256
`50d48b85f02a14f50614f0eea5cedbd7ccd68e9b2f7e1564fe8b41fcd2af9c4b`,
especially Sections 2, 5 and 7. This note does not use a new target event,
signed empirical law, ideal source law or post-cap independence.

## 1. Exact assumptions and cap-complement geometry

Fix n, x>0, t=n/x, a target w with ||w||_infinity=r in the retained
R073 target set, and any q in its COMPLETE matching interval. Use

    rho=tq, R=(3/2)q-h, a0=(1-epsilon)q,
    g=R-a0, kappa rho=rho+a0,
    ell_R(z_i)=min(1,(R+kappa rho-|z_i|)/(2R)),
    G(z)=product_i ell_R(z_i),
    bsoft(z)=min(1,theta/G(z)), theta>0.

Assume the existing active guards, in particular

    1<r<3/2, rho>4r, 0<q<=qmax<1,
    a0>0, 0<g<=R/3,
    0<J_i<=1+eJ,
    Nsum(q)/(2n ellpin(q))>=1-eN,
    Nsum=n+#{i: |w_i|<a0}, 2ellpin=1+a0/R.              (1)

Here J_i are the actual captured inverse-coordinate integrals from the
source proof. At this fixed q the full inverse-capture coefficient U is
the sum over EVERY captured signed source face, and the original proof
establishes U<=Uhi*n on the whole matching interval.

Inside the full source band q<1 and rho>1, m(z)=rho and
beta*rho=rho+q. Thus the original corner-cap complement is exactly

    D^c={z: z_i>=1-beta*rho=-rho+(1-q), every i}.         (2)

The source face with a NEGATIVE pin z_j=-rho has zero intersection with
D^c, since 1-q>0. A positive source pin z_j=rho satisfies (2) in its
pinned coordinate. The remaining conditions in (2) concern free
coordinates and are coordinatewise events.

## 2. Existing whole-matching count gives many fully captured lower strips

Write I(q)={i: |w_i|<a0} and mcount=|I(q)|. From g<=R/3,
a0/R>=2/3. Apply (1), retaining its exact calibration factor:

    n+mcount >= (1-eN)n(1+a0/R)>=(5/3)(1-eN)n,
    mcount >= nu*n, nu=2/3-(5/3)eN.                    (3)

This uses the already proved count bridge over the ENTIRE matching
interval. No q-grid or baseline-only replacement is made. It does not
assume an empirical law for signed target coordinates.

For i in I(q), w_i<=|w_i|<a0. Let d=1-q>0. Then the entire lower
strip B=[-rho,-rho+d) is captured by this coordinate. Indeed its lower
capture inequality is z_i>w_i-rho-a0, since w_i-a0<0. Its upper
capture inequality follows from

    z_i<-rho+d < w_i+rho+a0,

because d<1, w_i>=-r, rho>4r and r>1 imply
2rho+w_i+a0>2rho-r>7r>1>d. Boundaries have zero source-coordinate
Lebesgue mass, so strict and non-strict cap endpoints cause no loss.

## 3. A genuine pre-cap inverse-law estimate on each captured source face

Fix ANY captured positive signed source face j. Before cap restriction
and soft clipping, the source proof's normalized inverse-capture law
has independent free coordinates, with density

    1_capture(z_i)/(2rho*J_i*ell_R(z_i)), i!=j.

For each i in I(q)\{j}, the probability of its lower strip therefore is

    p_i = [1/(2rho J_i)] integral_B ell_R(z_i)^(-1) dz_i
        >= d/[2rho(1+eJ)],                              (4)

using ell_R<=1, its full capture just proved, and J_i<=1+eJ.
All integrals are actual finite-n coordinate integrals. At most ONE
counted coordinate is omitted by the source pin, regardless of which
target coordinate is pinned. If nu*n>=1, independence gives

    P_inv(D^c | face j)
       <= product_(i in I\{j})(1-p_i)
       <= exp[-(mcount-1)d/(2rho(1+eJ))]
       <= exp[-((nu-1/n)/(1+eJ))*((1-q)/(2q))*x].        (5)

The first inequality ignores all other free-coordinate cap conditions,
so it is an upper bound in the required direction. The negative signed
source faces have cap-complement probability zero. Set

    eta=(1-qmax)/(2qmax)>0,
    ccap(n)=(nu-1/n)*eta/(1+eJ).

Since q<=qmax, (5) implies on EVERY captured signed source face

    P_inv(D^c | j)<=exp[-ccap(n)*x].                    (6)

This is asserted only when nu-1/n>=0. If this positivity is not
verified, use the trivial probability bound 1 instead. Face coefficients
need not agree: multiply their nonnegative coefficients by the common
bound (6), then sum. One does not identify within-face and mixture
variance, nor assign a distribution to actual uncaptured increments.

## 4. Pointwise soft inequality cancels theta before the cap payment

For ALL source points, not only the retained event G>=theta,

    bsoft=min(1,theta/G)<=theta/G.                     (7)

At fixed q, the weighted captured cap-complement shell coefficient
is consequently at most

    theta * U(w,q) * exp[-ccap(n)*x]
       <= theta*Uhi*n*exp[-ccap(n)*x].                  (8)

This follows by applying (7), expressing the remaining inverse integral
on each captured source face using its pre-cap inverse law, applying
(6), and summing all faces. There is no assertion that independent
coordinates survive conditioning on D or soft weighting.

Let Lcap(w) be the loss from the full source row when imposing the
original source cap D, including both arad(w)=r/rhi and bsoft. The
matched source interval is q_+=(r+h)/(3/2), q_-=(r-h)/(3/2).
Integrating (8) against dq/q over THIS interval, not the entire source
band, gives

    Lcap(w) <= arad*theta*Uhi*n*log((r+h)/(r-h))
                         *exp[-ccap(n)*x].              (9)

Invoke the existing radial-width upper gate

    arad*log((r+h)/(r-h))<=2(1+eRad)h/rhi,
    k=n*h*theta/rhi, Mdegree=2(1+eRad)Uhi.

Then the complete local deletion bound is

    Lcap(w)/k <= Mdegree*exp[-ccap(n)*x].               (10)

It is uniform in every retained target w and the full matching interval.
There is no factor C_log/h and no exponent subtraction by chi. The
inequality (7) cancels theta against the theta in k exactly. It is
valid for every theta>0, including theta=exp(-chi*x).

## 5. Finite arithmetic illustration for the old R073 tuple at L>=4096

Keep its source band qmax=11569/12800, eta=1231/23138,
eN=2^-19, eJ=eRad=2^-20 and its Uhi. Assume the other local guards
in (1) and the radial-width gate have independently been proved at the
new domain. Use n>=10^6. Exact Fraction checks give

    nu-1/10^6 > 333/500,
    Mdegree < 3.

For x=(231/256)L, L>=4096, hence x>=3696, put

    cfixed=(333/500)*eta/(1+eJ)>0.

The new exact cap-only certificate verifies

    cfixed*3696>100,
    exp(100)>100^4/4!>3*2^20.

Equations (10) and these inequalities imply over this whole local tail

    Lcap(w)/k < 3*exp(-100) < 2^-20=eS.                (11)

The simple base n>=10^6 follows, for example, from the positive-series
lower bound n=e^L>=L^64/64! and its exact value at 4096. All displayed
finite gates are checked by the new verify.py. This is not a rerun or
revision of the R079 batch's 53-gate receipt, and does not certify its
other assumptions or its tau payment. This illustration is deliberately
coarse: the cap exponent lower envelope at 4096 is about 130.96, but
the proof needs only the rational threshold 100.

## 6. What changes, what remains a bottleneck, and acceptance scope

The old deletion divided its unweighted mass by k, costing a band-width
factor of order L and requiring the competing exponent eta-chi>0.
The new deletion estimates the weighted capture integral directly and
uses the existing count event to extract a product small probability.
It removes that particular eta-chi restriction. A useful positive eta
is still needed: taking qmax all the way to 1 makes ccap vanish, so
this lemma does not make arbitrary upper source bands free.

The price is exactly the already available count gate. Without a
lower bound on mcount, a universal exponential decay cannot be claimed:
the inverse-capture law can remove lower strips entirely in coordinates
with w_i>=a0+d, making their cap-complement probabilities one. Existing
source/target profile conditions must not be silently replaced by a
uniform signed-coordinate law. Here (3) supplies the needed lower
count and the arbitrary omitted-coordinate cost is paid by the -1.

This cap improvement does NOT reduce tau/k. Its rate
sigma/2-chi(1-sigma), the Bernstein retention gap chi-m, target mass,
row/column constants, finite-a same-history transfer and complete root
cost all remain separate gates. One must not activate a new operator
bound solely because (11) succeeds.

The final kernel, both weights <=1, original source cap and actual
history are unchanged. Thus, once independently accepted, (10) may
replace only the source-cap payment in a full theorem; its finite-a
and arbitrary-L1 bridges require the same existing proof review. It
does not invoke an ambient probability law for actual increments.

The independent sparse clone from R078 is reused. Only
`upper/round079/weighted_cap/` was written, with no commit or push.
The source hash is checked before the new arithmetic. Python was
actually run for the final certificate; no shared MATLAB, parameter
grid, CP-tail calculation, lower-bound work or n=7 work was used.
