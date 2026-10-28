# R073 full-chain independent audit — ACCEPT

Author source, read in full:

`/Users/zhengzhihao/Desktop/T/cube-r040-cover-selector/cube-dim-7-bounds/upper/round073/cantelli_frontier/PROOF.md`

SHA256 `50d48b85f02a14f50614f0eea5cedbd7ccd68e9b2f7e1564fe8b41fcd2af9c4b`.

The complete mathematical chain is accepted for every integer n with
L=log n>=8192:

    C_n <= F(A_n(t), B_n(t)-(291/1000)n)
        < n log n-n log log n-(47/250)n,
    t=n/[(231/256)L].

This is the original centered Lebesgue cube maximal operator, arbitrary
L1 input, all positive witness radii and the same genuine finite-history
algorithm. It is not merely a scalar/frontier model. The earlier model
threshold 40448/5, a grid-optimum identity, or a global optimality claim
is not part of this acceptance. The leading coefficient remains one.

Only the owned private audit directory was written. No author/main
source, frozen R071 material, Git state or MATLAB session was modified.
No production, grid, n=7 or lower-bound computation was performed.

## 1. New tuple and local domains

The source correctly distinguishes qmax=11569/12800 from
1-sigma=231/256. The former controls the source band and eta; the
latter controls x, A and the scalar price. Their resulting cap and
tau decay rates are different and were not inherited from R071.

At 8192, delta/h=12; thereafter h decreases. The full matching interval
and the baseline q_b=q_+ stay in the source band. Active R lies in
[r-2h,r], exceeds one, and r<3/2. The new finite guards give
epsilon q<=h/s, rho>4rhi, kappa>1 and beta<2. Thus 0<g<=R/3.
Unmatched sources with R<1 have zero final column. The local bounds
are applied only on their genuine relative-gap domain, not by renaming
the earlier absolute-gap domain or by importing R070's different
h-coefficient interval.

The frozen local statements were checked against their actual guards:
R064 c60d48d1..., R066 67bdf01d..., and R070 4d2f4fc2.... In particular
the strong constants 3/rho, 1/(2rho), 1/(10rho) really occur in those
proofs. They are integrated a.e. using AC across positive-part kinks.
The correct mean-one baseline is R=r at q_b, not q=r/s.

## 2. Target laws, pins and six bad events

The exact J,J1,J2 tail integrals in section 3 retain all prefactors and
central deletion terms. J2 is the second integrand profile, not J1^2
or the second moment of the whole sum. Target Fubini cancels capture
probability ell_R against one inverse ell_R and yields full-source
phi and phi^2 means. It does not equate target and inverse-source laws.

The lower finite-h J2 mean genuinely uses
a-5/6 in [0,h/(3R)] and R/q=s-h/q. Its lower bound would not follow
from g/R<=1/3 alone. The source pays n-1 free coordinates and the
deterministic target pin separately. The exact upper pin
1/(250rho_b) is used only for a pinned target, not an arbitrary source
pin. New independent alternating log8/9 enclosures reproduce the
claimed [.00528,.005295]x mean corridor for the actual new h.

The six target conditions are imposed on each full signed ambient
face. Free coordinates are independent there; the six events need not
be mutually independent. Baseline log mean/variance, both counts,
target cap, J1 and the two-sided J2 condition each cost <=epsA/8.
For J2 both excluded tails lie inside one Chebyshev event, so there is
no missing factor two. For J1 the deterministic pin is at most half
the 10^-5 x slack at n>=9000, leaving the stated Chebyshev gap.

The cap-complement rate (3/4)^(n-1) is per positive pinned face;
negative faces already lie in D. No erroneous one-half normalization
is used. The union leaves >1-epsA of each shell, giving target mass
>(3043/10000)n. The exact logarithmic mass check was recalculated in
the new preflight, rather than copied from the old tuple.

## 3. Whole-q control and facewise Cantelli

The sign ellpin'<0 and monotonicity of Nsum give the correct endpoint
ratio calibration. The new 2h+1/n tolerance is paid, and hence the
ratio stays in 1+-eN on the entire interval. The log-product shift is
6Lh=6/768. It is explicitly deducted in P0 using exp(-u)>=1-u.
The individual J_i denominator bounds are pointwise and are separate
from this n-fold product bound. Summing every available source sign
then gives exactly U in equation (17), with Ulo/Uhi.

Strong derivative integration transports J1 by <=xh and J2 by <=xh/5.
The new h pays both margins. On each fixed captured signed source
face, the pin in V=-log G is deterministic and <=.183; free inverse
coordinates have exact means J1_i/J_i and raw second moments J2_i/J_i.
Subtracting their squared means gives the true conditional variance.
For its upper bound, adding back the omitted nonnegative J2 profile
and dividing by 1-eJ gives <=.00531x. No pin variance, lower variance,
normal approximation, mixture variance or actual-mu randomness is
asserted.

The exact base retention exceeds 9717/10000. With a=chi-m_*>0 and
u=a x-p>0, the derivative of u^2/x is a^2-p^2/x^2>0. Thus the fixed
retention is uniform over every larger L, not just a rounded base
number. Cantelli is applied separately in each source face, then
the probabilities are mixed. The bsoft/(theta/G) ratio proves the
soft-row retention. All this occurs before cap conditioning.

## 4. Complete radial degrees and tight cap normalization

The full-band bsoft is defined before D removal; arad and bsoft are
both bounded by one. Radial cancellation is exact, with only the
atanh(h/r)/(h/r) width correction. The normalizer retains rhi:
k=nh theta/rhi. Raw row bounds therefore have precisely the stated
2Ret(1-eRad)Ulo and M coefficients.

For the column, all target signs are retained. Since kappa rho>|z_i|,
at least one sign is available, and T_i/(2ell_i(r))<=1. On matching
r>=R, so product ell_i(r)<=G(z). Integrating the ordinary matched
width <=2h after radial cancellation yields column <=2k; inactive
columns are zero. M>2 consequently bounds both degrees.

The cap loss is paid by its UNWEIGHTED ambient mass only after these
full-band laws have been established. Its logarithmic band mass is
n log(qmax/qmin); the older linear upper bound and factor two are not
silently substituted. eta=1231/23138 gives the correct exponent.
The finite-n exponential correction is positive and explicitly paid:

    eta*x/n<=1/1025,
    exp(eta*x/n)<=1/(1-eta*x/n)<=1025/1024.

Dividing by k preserves C_log*rhi/h. The exact rates displayed in the
proof are distinct and positive:

    cS=13627383/3790929920,
    cT=1447/327680.

Independent eight-term atanh-log bounds and P40 already verified both
base payments, with the same correction factor and complete prefactors.
No rounded exponential or logarithm is needed.

## 5. Entire L tail and the complete root buffer

For a positive prefix P_m of degree at least two, omitted terms in
yP_m'-P_m and yP_m'-(3/2)P_m are nonnegative, giving the respective
lower bounds -1+y^2/2 and -3/2-y/2+y^2/4. The two new rates separately
have cL0>4. Thus the cap L/P_m(cS L) and tau L^(3/2)/P_m(cT L)
bounds decrease throughout the tail. The rational 91 bounds sqrt L
only at L0; it is not used as a uniform bound on sqrt L.

All remaining adverse envelopes decrease after substituting N(L).
In particular x/N, epsilon q/h, probability rates, individual-J and
count/radial errors improve; source lower radii increase. The finite-h
J2 lower mean increases, its upper decreases, and the log-product loss
is constant. Retention has its own explicit monotonicity proof. There
is therefore no residual unquantified eventual domain in the result.

The genuine B lower bound B>=nL/4 makes B-gamma n and tau positive.
The exact log-price bound uses sigma, not qmax, and is <257/2500.
The identity gamma-price-47/250=1/5000 leaves a positive budget for
all A/cross/self terms. A new P16 check verifies the full remainder
1/n+exp(-sigma L)+2sqrt L exp(-sigma L/2)<1/5000 at L0. The same
prefix inequalities cover powers zero and one-half on the entire
tail. Thus the strict full F inequality follows, not merely a
comparison of gamma coefficients.

## 6. Genuine finite history and general L1 transfer

A and P are point-disjoint. An actual target item meeting A has
r_j<3, whereas a source point belongs to an item with r_i>=rho>4.
Hence each included source precedes the same target item in the
actual algorithm. With algorithmic a>=8n^4, x>1 gives t<n, beta<2,
c_cap>1 and zeta-(beta+2)/(a+1)>0. Included pairs really capture that
target center in its own eligibility test. Included increment pieces
are paid by original item weights, without copying a whole item into
every band. The bounded weights leave actual row <=arad*tau<=tau.

The finite envelope equals c_a lambda on A and P, with c_a=alpha^-n>=1,
and D lies in the finite cap. Energy, degrees and target mass scale
by c_a^2,c_a,c_a respectively. The two deficit terms are paid by
nu_a(A)+nu_a(P)<=Delta_a because supports are disjoint. This gives
the displayed Gamma n deficit uniformly over every genuine finite
history. No ambient coordinate law is assigned to mu_a.

The full cover factor, terminal self term and scalar a-limit are
retained. At tau=sqrt(A(B-gamma n)), their limit is exactly F(A,B-gamma n).
Positive witnesses, compact finite subcovers and inner regularity
provide arbitrary L1 and all positive terminal radii, with no compact
support, common selected history across unrelated inputs, selector
convergence or scale truncation. Thus the actual operator conclusion,
not merely a relaxation, follows.

## New computation and stage identities

* Preflight stage: seven independent Fraction groups, actual PID 1890,
  PASS. N32, log8 plus geometric tail and P40. It covered sensitive
  retention, mass, gamma, cap/tau and price margins before source freeze.
* Frozen-source stage: six additional groups, actual PID 2503, PASS.
  N32, alternating log8/9 and P16. It covered the previously unchecked
  local geometry, target probability budgets, actual finite-h J2 mean,
  profile transport and complete positive root buffers.

The two stages and their receipts remain distinct. Neither was rerun.
The author 53-gate program and R071 programs were not run. No new
search was performed. Finite checks support the exact margins; the
analytic arguments above establish all-tail and measure-law transfer.

Final source identities and audit-file hashes are in sourcehash.json.
