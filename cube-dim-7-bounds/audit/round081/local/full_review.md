# R081 full s=4 eventual theorem — independent ACCEPT

The author's frozen PROOF.md, verify_constants.py and verification.json
were read in full. This review accepts the complete eventual statement

    there is a finite L_* such that, for every integer n with
    L=log n>=L_*,
    C_n<=Phi81=F(A_n(t),B_n(t)-(13/10)n)
       <T-(3/8)n,
    t=n/[(2/5)L], T=n log n-n log log n.

This is the original centered axis-parallel Lebesgue cube maximal
operator for arbitrary L1 inputs, every positive strict level and every
positive witness radius. **No explicit L_* is certified. In particular
this is not a theorem at L=4096.** The exact comparison to Phi80 holds
on the actual common domain, L>=max(L_*,4096).

All source files were read only. This staging directory is not a clone,
and new files stay in the owned `audit/round081/local` directory.
No original checkout, author file, Git state, shared MATLAB, paper or
accepted ledger was changed. No commit or push was made. The author's
13-check program was not run.

## 1. New local lemma, rederived without s=3/2 constants

`proof.md` in this directory independently derives all three EXACT
coordinate formulas, the value bounds and the strong AC derivative
estimates from their one-dimensional integrals. It explicitly allows
the central-deletion region |w|>R. It supplies the sufficient local
domain h<=1/1000, rho>=12, 1<=R<=4 and the stated relative gap.
These concrete additions hold eventually in the proposed geometry.

The independent bounds are stronger than those used by the author:

    |J-1|<=4/rho, |d_q log J|<=711/(20rho)<100/rho,
    J1<=1/rho, |d_q J1|<=5601/(1000rho)<40/rho,
    J2<=1/(3rho), also |d_q J2|<3/rho.

For the fixed target, t,h,epsilon are held constant during q
differentiation. The positive parts are Lipschitz, Q,F are smooth on
the compact relative domain and J is bounded below by 2/3. Thus the
profiles and log J are AC, and a.e. estimates may be integrated across
all onsets. The proof does not rely on a finite difference or grid.

For |w|>R the deleted coordinate interval has length m>g. It removes
the near edge plus a central interval, where phi=0 and ell=1, without
reaching the opposite edge. Therefore J retains the necessary -m
term and J1,J2 simply have their near-edge v contribution set to zero.
This resolves the potentially problematic full-matching region.

The bounded local Python/SciPy run PID 28899 is supplementary evidence:
24 fixed points spanning full capture, edge deletion and central
deletion reproduce all three integral formulas to 1.13e-16 absolute
quadrature error. Exact rational inequalities support the constant
bookkeeping. These computations do not establish the uniform theorem;
the analytic local proof does.

## 2. Actual matching domain and five target events

The new fixed tuple is s=4, q in [1/4,.99], sigma=.6, chi=.7,
h=L^-2, delta=10^-5 and all errors e=epsA=10^-4. Source radius
rho=tq grows uniformly, while h and epsilon*q/h tend to zero.
Thus the full matching q_±=(r±h)/4 lies in the full source band once
2h<delta. Its active R=4q-h lies in [r-2h,r] and exceeds one.

The exact relation

    g-3R/4=epsilon*q-h/4<=0

pays the new relative-gap domain uniformly. On the full band R>0,
a0>0 and G>0 eventually; no active R>=1 estimate is applied to
unmatched edge sources. Both bands lie within norm<t, and rho>4rhi
separates them. All these statements follow from displayed explicit
envelopes, not an inherited old threshold.

At q_b, R=r. Averaging capture over the full uniform target coordinate
gives exactly ell_R(z), hence E_W J=1. If |J-1|<=6/rho<=1/2,
the elementary inequalities log(1+y)>=y-y^2 and
|log(1+y)|<=2|y| give

    -36/rho^2<=E_W log J<=0,
    Var_W(log J)<=144/rho^2.

The target pin is deterministic and has |log J|<=12/rho. Thus the
full-face baseline sum has the author's mean and variance envelopes,
uniformly tending to zero because n/rho^2<=16x^2/n and 1/rho<=4x/n.
Chebyshev supplies its e tolerance eventually.

Full-target Fubini for J1 gives the exact full-source loss mean

    E_W J1=(2R/rho)[1-a+a log a], a>=5/8.

The bracket decreases as a increases, and R/q<=4. Thus the free
sum mean is <=mu4*x, mu4=3-5log(8/5)<.65. The actual target pin
adds at most 1/rho, eventually <=.005x, giving mean <=.655x.
The free-sum variance is <=n/rho^2<=16x^2/n. Chebyshev across the
remaining .005x gap gives exactly the stated 640000/n bound for
sum J1>.66x. No second-profile target event is needed.

At each endpoint the pin's excess count indicator is zero because
a0<r. The other n-1 indicators are independent Bernoulli(a0/r), so
each count failure is <=4/(e^2 n). For the original target cap,
negative pinned faces lie inside D once t>rhi; positive pinned faces
have complement fraction 1+1/(2t)-1/(2r)<.9 eventually since rhi<4.
Their cap failure is <=.9^(n-1).

These are exactly five target events: one log, two counts, D and J1.
Each failure is eventually <=epsA/8, uniformly over every radius and
signed target face. They need not be independent of one another.
The union leaves mass fraction >=1-epsA and gives

    lambda(A)>=(1-epsA)n log(rhi/rlo)>1.37n.

All probability laws are ambient full-target laws before conditioning,
not assumptions on actual increments.

## 3. Whole matching, source mean, TRUE variance and eventual retention

The interval length is h/2. The new log derivative bound gives

    sum_i |d_q log J_i|<=400x,
    whole-interval product-log shift<=200xh=80/L.

For L>=800000 this shift is <=e. Together with the baseline log
event the product lies in [P0,1/P0], P0=1-2e, using
exp(-2e)>=1-2e and exp(2e)<=1/(1-2e). Individual J_i lie in
[1-e,1+e] eventually, separately from the n-fold product estimate.

The derivative ellpin'=-(1-epsilon)h/(2R^2)<0 and the nondecreasing
count give endpoint calibration of Nsum/(2n ellpin). The q_+ mean
ratio is <=1. At q_- its shortfall is <=2h+1/n, using a0<=1 and
R>=1. Count fluctuation contributes <=e/4. Once 2h+1/n<=3e/4,
the ratio lies in [1-e,1+e] over the WHOLE interval. Summing every
captured signed source face and omitted coordinate yields the stated
Ulo*n<=U<=Uhi*n. There is no baseline-only substitution.

The new J1 derivative bound integrates to at most 80xh. Once
80h<=.01, the whole-matching first profile is <=.67x. In each
captured SOURCE face, before cap deletion or soft clipping, the free
inverse coordinates really form a product law with first moments
J1_i/J_i and raw second moments J2_i/J_i. The source pin is
deterministic and <=log(8/5)<.5. Thus

    E V<=.5+[.67/(1-e)]x<=.5+.671x,
    Var V=sum_(i!=j)[J2_i/J_i-(J1_i/J_i)^2]
       <=n/[3rho(1-e)]<=4x/[3(1-e)]<2x.

This is the TRUE centered variance. The pin contributes no variance;
J2 is neither J1 squared nor a second moment of the whole face sum.
The denominator and the arbitrary omitted coordinate are handled in
the correct direction. A second-profile good-target event is not
required for this coarse pointwise variance bound.

For u=.029x-.5>0, the bad event V>.7x implies V-EV>u, so
Chebyshev gives probability <=2x/u^2, uniformly in target, matching
q, face and sign. This tends to zero as x grows and is eventually
<=.001. Equivalently u^2/x grows without bound; no numerical normal
approximation, lower variance or Bernstein hypothesis is imported.
After the uniform facewise estimate is established, faces may be
mixed. The soft ratio min(1,G/theta) equals one on V<=chi*x and is
nonnegative elsewhere, giving the claimed raw Ret=.999 lower row.

## 4. General-s weighted cap and complete radial degrees

Since a0/R>=1/4, the whole-matching count yields

    mcount>=nu*n, nu=1/4-(5/4)e=.249875>0.

On q<1 the original cap complement is exactly all
z_i>=-rho+(1-q). Negative source pins have zero cap-complement mass.
For any positive pinned source face, each counted FREE coordinate
captures its entire lower strip, because w_i<a0 gives the lower
endpoint and rho>4rhi gives the other endpoint. Its inverse-law
probability is at least (1-q)/[2rho(1+e)]. An arbitrary source pin
omits at most ONE counted coordinate. Product-law independence gives

    P_inv(D^c | face)<=exp[-(nu-1/n)*eta*x/(1+e)],
    eta=(1-qmax)/(2qmax)=1/198.

The positivity nu-1/n>=nu/2 holds eventually. Pointwise
bsoft<=theta/G converts this to a bound for the weighted captured
loss; the estimate is taken BEFORE cap conditioning. Summing all
source-face coefficients and integrating only over the actual matched
dq/q interval gives the same radial log-width as before. Division by
k=n*h*theta/rhi cancels theta exactly:

    caploss/k<=M exp[-(nu-1/n)*eta*x/(1+e)] ->0.

This is the new weighted deletion, with no unweighted source-band
mass, no chi subtraction in its exponent and no obsolete eta>chi
restriction. The loss is eventually <=eS=e. Both kernel weights
remain <=1; the probability estimate is never imposed on mu_a.

The exact raw radial row normalization is unchanged by s:
arad*log((r+h)/(r-h)) is within 1±e of 2h/rhi eventually. It gives
raw row >=2Ret(1-e)Ulo*k and row <=M*k.

For the complete column and matched r>=R, kappa*rho-|z_i|>=a0>0
ensures ell_i(r)>1/2. The expression
min(1,1/2+(kappa*rho-|z_i|)/(2r)) decreases in r, so
product ell_i(r)<=G. Every target pin sign is counted and
T_i/(2ell_i(r))<=1. Thus the complete target-face coefficient is
<=nG. The arad factor cancels dr/r, the ordinary matched target
width is <=2h and bsoft*G<=theta. The entire column is <=2k<M*k,
with inactive columns zero. The enlarged target radius introduces no
extra degree factor.

## 5. One finite eventual threshold, tau and gamma greater than one

All new geometry, local and five-event envelopes are uniform in the
full stated target/source domains. The error budgets are fixed and
positive. Their adverse polynomial-over-n bounds tend to zero,
80/L and h tend to zero, the source cap bound decays exponentially,
and the Chebyshev retention bound tends to zero. There are finitely
many such inequalities, so their finite eventual thresholds have a
finite maximum. This is a rigorous existence argument; it does not
compute that maximum or claim it is 4096. The local h<=.001 and
rho>=12 conditions from this independent proof are among those
eventual guards.

The original cap identity gives B_n>=n log(t)/2. Since log t>=L/2
eventually, B_n>=nL/4>1.3n. Therefore tau=sqrt[A(B-1.3n)] is
positive. The restriction gamma<1 sometimes appearing in old finite
gate programs is not needed by the operator argument. Here it is
replaced by this direct positive-root proof.

A<=exp x and B<=nL give

    tau/k<=rhi L^(5/2) exp(-L/50) ->0,

with the strictly positive rate sigma/2-chi(1-sigma)=.02. The factor
L^(5/2) correctly includes h^-1=L^2 and sqrt L. This same numeric
tau is used throughout the actual selection history, not a separate
post-selected threshold. Paying this tau and weighted cap leaves
row-tau>=b*k, while both degrees remain <=M*k.

The fixed rational deficit inequality is (1-e)Gamma>1.3, with
Gamma=mL*b/M about 1.367125. Gamma>1 is consistent: the total full
corner-envelope mass is B of order n log n and there is no universal
unit-n upper limit on the local missing mass in the energy argument.

## 6. Original finite-a history, arbitrary L1 and all positive scales

For every a_alg>=8n^4, t<n eventually implies a_alg>2t and ccap>1.
An actual target item meeting A has ccap*r_j<=||u-w||<alpha+rhi<5,
so r_j<5. An actual source increment has r_i>=rho>4rhi>5. Hence
it precedes that SAME target item in the original largest-radius
history. The new target bound is 5, rather than the old 3, and the
source bound explicitly exceeds it; no old separation constant is
silently reused.

The displacement guard

    alpha*beta-kappa-2/(a_alg+1)
      =n^-3-(beta+2)/(a_alg+1)>0

holds with beta<2. Since r_j<r_i, included true centers lie within
[kappa+2/(a_alg+1)]r_i<alpha*beta*r_i. They really capture that
target in its own original eligibility test. Integrating only actual
included increment pieces costs at most their genuine item weights.
The two bounded kernel weights give actual row<=arad*tau<=tau.

The finite envelope equals c_a*lambda exactly on A and P, because
D lies inside its cap by the displayed positive algebraic margin.
Energy, degree and target mass scale by c_a^2,c_a,c_a. The two
energy-difference missing-mass terms cost at most
c_a*M*k*[nu_a(A)+nu_a(P)]<=c_a*M*k*Delta_a, since A and P are
disjoint. Actual energy is <=tau*c_a*lambda(A). Therefore

    Delta_a >=[c_a^2 E-tau*c_a*lambda(A)]/(c_a*M*k)
       >=[E-tau*lambda(A)]/(M*k)
       >=lambda(A)*b/M>=Gamma*n>1.3n.

The use of c_a>=1 has the right direction. Every step remains valid
for gamma>1; none divides by 1-gamma or uses a deficit capped by n.
This is uniform over every original finite history and finite dilation.

For general L1 use |f| and true positive-radius witnesses for every
point of a strict level set. Open shrunk witnesses have finite
subcovers of compact subsets; selected dilates remain true witnesses.
The terminal-radius normalization permits every positive radius.
Lebesgue inner regularity and level homogeneity recover the complete
operator statement without compact support, cutoffs, scale lattices,
actual-increment independence or selector convergence.

## 7. Full root cost and complete comparison with Phi80

The ENTIRE finite cover factor and terminal self-weight are retained.
Taking only the numeric a_alg limit, already uniform over witness
families, yields exactly

    C_n<=(1+A/tau)(tau+B-1.3n)=Phi81.

The original scalar bound B<=T+n log(5/2)+1 and
A+2sqrt(AB)<=n exp(-.6L)+2n sqrt L exp(-.3L) then give the complete
strict display. The sum of the positive A/cross/self remainders tends
to zero and is eventually <1/200. Since

    log(5/2)<23/25,
    13/10-23/25-3/8=1/200,

one obtains Phi81<T-(3/8)n. No positive root term has been omitted.

For the old R080 function, the exact original integral bound
I_n(t)<=2t/(n-1), its sigma=.11 and gamma=.39029 give

    Phi80>=B80-.39029n
       >=T-n[.28029+4/L].

This is a lower bound for the old upper-bound FORMULA, not for C_n.
It remains valid on R080's accepted positive-root domain L>=4096.
The strict separator 3/8>.28029+4/4096 therefore proves
Phi81<Phi80 whenever both branches are valid, namely
L>=max(L_*,4096). The comparison does not create an explicit threshold
for R081, extend any old branch, or claim global optimality.

## 8. Acceptance object and actual evidence

This full review, the independent coordinate proof, and the unchanged
author proof together close the local, probabilistic, all-tail,
finite-history and operator interfaces. No mathematical gap was found.
The author's receipt is only a fixed-constant calculation; it is not
being relabeled as a finite-threshold computation.

`full_check.py` independently checks sensitive fixed gaps using exact
Fraction arithmetic and log16 enclosures. It does not import or run
the author program and does not certify an explicit L_*. The local
quadrature/finite differences remain diagnostic evidence only. Source
hashes and the actual new-run identities are recorded in the receipts
and manifest; the original source hashes are rechecked on completion.
