# R051 flexible-budget independent audit

Conclusion: **ACCEPT**, for each complete fixed tuple in the stated domain, with a tuple-dependent sufficiently-large-n threshold. No explicit finite threshold, n-dependent adaptive parameter choice, continuous global optimum, or numerical-grid proof is accepted by this conclusion.

Audited frozen source:

    /Users/zhengzhihao/Desktop/T/cube-r030-geometry/cube-dim-7-bounds/upper/round051/flexible_budget/PROOF.md
    SHA256 aeddd67b7a2109d402ea76d9f6fb3d93accf44bea80240340bd163759a7c9df6

The author supplied this exact source as ready and frozen. This auditor did not edit it. Work here is confined to the private clone's `audit/round051/flexible_budget/` directory. No production search, old grid rerun, n=7 work, MATLAB or Git write operation was used.

## 1. Accepted parameter and budget statement

The full domain is necessary:

    a>3, s>2a, 0<d<=1/16, d<1/(32a), sigma=sd^2<1/2;
    0<epsJ,epsA,epsD<1, epsS,epsT>0;
    bracket=(1-3/a)*Ulo-epsS-epsT>0.

With P0=3579/4096,

    Ulo=P0/[2(1+epsJ)],
    Uhi=1/[P0(1-d/2)(1-epsJ)],
    mass=(1-epsA)*d/[2(1+3d/4)],
    M=max(4,2Uhi/(1-3d/4)),
    Gamma=mass*bracket/M,
    gamma=(1-epsD)*Gamma,
    c_raw=gamma-sigma/(1-sigma).

The audited theorem proves eventually Cn<=F(A0,B0-gamma*n). Every **fixed** 0<c<c_raw then gives Cn<T-c*n, where T=n log n-n log log n. The endpoint c=c_raw is not included by the displayed o(n) envelope.

The proof does not use a sign-reversing shortcut. Positive bracket, positive mass and M>=4 imply Gamma>0. The domain also gives Ulo<1/2, mass<d/2 and bracket<Ulo, hence Gamma<d/16<=1/256. This justifies the quoted small deficit independently of any optimized tuple.

## 2. Coordinate proof survives the relaxed d domain

The old tiny fixed-d restriction was not needed for the local inequalities. On the full baseline-to-match q interval, eventually R,q>=1/2, g<=d<=1/16, T<=2d, rho>=4R and g'=zeta*t<=d. At baseline R=r>1. The baseline q may be h outside the original source-band endpoint; it still lies in the same fixed near-one interval and satisfies every local guard.

For f(y)=2R/(2R-g+y), 0<=y<=g, direct differentiation gives

    |f_R|=2(g-y)/(2R-g+y)^2<=8d,
    |f_g|=2R/(2R-g+y)^2<=4,
    f<=2.

These bounds already hold at R=1/2, g=d=1/16, and improve as R increases or g decreases. Integrating the appropriate branch of A and differentiating its endpoint gives the source's conservative bounds

    |A_R|<=16d^2, |A_g|<=8d, A_T<=2,
    |B_R|<=8d^2, |B_g|<=4d.

For T>=g the endpoint derivative cancels because f(g)=1; no derivative jump was omitted there. Capture onset T=0 is locally Lipschitz, so an a.e. bound may be integrated.

In differentiating J=1+(2B-A)/(2rho), the numerator derivative contributes at most `(24d^2+I_tail)/rho` using the displayed coarse constants. The derivative of 1/rho must also be counted. Since q>=1/2, B<=g^2/(3R), A=0 off the tail and A<=4d on it, that extra term is bounded by `(4d^2/3+4d*I_tail)/rho`. With d<=1/16 their sum is strictly below `32(d^2+I_tail)/rho`. As J>=1/2 eventually,

    |partial_q log J|<=64(d^2+I_tail)/rho.

Here onset implies |w|>r-2d: along the interval |r-q-d|<=h and zeta*rho<=h/2, with h<d/8. Thus the tail indicator used in the source is valid. No omitted O(d) bulk term remains from differentiating rho.

For the baseline free coordinate, EJ=1 follows exactly by Fubini. Using B<=g^2/(3R), A(T)<=4T/3 and ET^2=g^3/(3R) gives the stated conservative `E(J-1)^2<=3g^3/rho^2`. The pointwise bound |J-1|<=g/rho<=1/8 then yields E log J>=-6g^3/rho^2 and Var log J<=12g^3/rho^2. The pinned target coordinate is deterministic and contributes O(g/rho); it is not assigned an independent free law.

## 3. Good-target mass and product budget are different limits

At a fixed target face/radius, the n-1 free coordinates are independent uniforms. Their logarithmic sum has absolute mean and variance O(nd^3/t^2)=O_d(L^2/n)->0 uniformly in the radial interval; the pinned term tends to zero. Thus the baseline |sum log J|<=1/8 event has probability 1-o(1).

Also Kw=1+Bin(n-1,2d/r). Its mean is at most 1+2nd and variance at most 2nd, so the probability of Kw>4nd tends to zero for every fixed d>0. The original target corner complement is exponentially small. These uniform conditional estimates justify the arbitrary fixed epsA loss and exact mass lower bound

    lambda(A)>=(1-epsA)*n*log[(1+3d/4)/(1+d/4)]>=mass*n.

For every target already in this good set, summing the log derivative, n/rho<=2L, Kw<=4nd and integrating the q interval of length at most 2h gives 1280Ldh=5/4096. Together with 1/8=512/4096 this gives |sum log J|<=517/4096. The exact inequalities exp(-z)>=1-z and exp(z)<=1/(1-z) yield

    3579/4096 <= product J <=4096/3579.

Separately, the **individual** pointwise J tends to one uniformly for all allowed coordinate values and q. This permits any fixed epsJ eventually. It is not used to raise `(1+-epsJ)` to the n-th power. The fixed product enclosure is retained independently, which avoids precisely that invalid shortcut.

## 4. Full inverse-face sum and correct soft-retention probability

The shell coefficient is the actual face sum

    U=(product J)/(2ell_*) * sum_i N_i/J_i.

Pinned source ell_* is fixed on each face, `1-d/2<=ell_*<=1`; 1<=Ni<=2. Combining the fixed product enclosure with the individual J enclosure gives exactly Ulo*n<=U<=Uhi*n. This is not a conclusion from a single deleted-coordinate product.

The soft weight b=min(1,theta/G), theta=exp(-aLd^2), satisfies 0<=b<=1 and b<=theta/G. The lower retention is correctly evaluated under **each captured face's inverse-row law**, whose free density is proportional to capture/ell and whose normalizer J is at least 1/2. Because ell>=3/4,

    E_inv[-log ell]<=2 E_full[(-log ell)/ell]
                       <=8g^2/(9Rrho)<=d^2/rho.

The pinned coordinate contributes at most d. Therefore every captured face, and their normalized inverse-row mixture, has E_inv[-log G]<=d+nd^2/rho<=3Ld^2. Markov gives tail probability at most 3/a. It is essential that this is the inverse-row law, not the raw source law.

On G>=theta, b=theta/G; on the tail it remains nonnegative. Thus the lower soft coefficient `(1-3/a)Ulo*n*theta` and upper `Uhi*n*theta` have the correct directions. Ordinary matching width lies in [h,2h], and source log-width lies in [h,2h/qmin], giving the exact row constants stated in the source.

## 5. Both exponential payments and degree budgets

The source cap complement has full unweighted mass O(nd*exp(-Ld/32)). Indeed on a positive pinned source face each free coordinate has capture fraction 1-v/(2rho), v>=d/4; the negative pinned face is excluded. Using n-1>=n/2 and x>=L/2 gives the stated exponent. Since b<=1, this total mass bounds the lost weighted row for each target. Its ratio to k=nh*exp(-aLd^2) is polynomial in L times

    exp[-(d/32-a*d^2)L].

The strict d<1/(32a) makes this tend to zero, hence any prescribed fixed epsS is available eventually. This uses actual unweighted cap mass; no favorable unknown-mu distribution is assumed.

The column estimate is pointwise: target lengths decrease with target radius r>=R; omitted-coordinate face products cost at most the displayed factor two; matching log-width is at most 2h/R. Thus Craw<=4nhG/R, and multiplying by b gives at most 4k. The row upper bound is 2Uhi*k/qmin. Their maximum is exactly bounded by M*k.

The selection threshold uses A0<=n*exp(-sigma L), B0<=nL, so

    tau/k <=(sqrt L/h)*exp[-(s/2-a)d^2 L]->0.

The strict s>2a is necessary for this displayed argument and is correctly present. There is no missing theta^-1 factor: b<=1 preserves the actual eligibility row directly. The source cap condition and the tau condition pay different exponents and cannot replace each other.

## 6. Actual mu, one history and the general operator

For any original finite history the target/source kernel uses the same actual uncaptured-increment measure mu_a and corner envelope lambda_a, with mu_a<=lambda_a. Target and source regions are disjoint. On both, lambda_a=c_a*lambda exactly, c_a=alpha^-n>=1, by the source's D-subset-H_a computation.

For every algorithmic dilation mathfrak a>=8n^4, the center-displacement guard is uniform eventually: kappa+2/(mathfrak a+1)<alpha*beta. Actual target items meeting A have radius below 3; actual source items have radius above 4, so each included source precedes **the same actual target**. The true center-capture argument then bounds its original eligibility row by tau. Multiplication by b<=1 preserves that bound, including all radial and target-good restrictions. No independent history or regularity of mu is introduced.

For nu=lambda_a-mu_a, the nonnegative kernel identity is

    E_lambda,lambda-E_mu,mu
      =integral K dnu_A dlambda_P+integral K dmu_A dnu_P
      <=D_a[nu(A)+nu(P)]<=D_a*Delta.

Together with E_mu,mu<=tau*lambda_a(A), the energy, degree and mass scaling by c_a^2,c_a,c_a gives a deficit at least

    lambda(A)*[(1-3/a)Ulo-epsS-epsT]/M>=Gamma*n.

Since gamma=(1-epsD)Gamma<Gamma, the chosen gamma*n is a uniform valid finite-history deficit for this fixed tuple. The actual cover/load inequality then holds for every sufficiently large mathfrak a, and only these scalar inequalities are passed to the mathfrak a->infinity limit.

The original finite compact-level-set cover is unchanged: every selected dilated cube remains a witness of the original arbitrary L1 input; there is no bounded-scale replacement. Lebesgue inner regularity supplies the whole level set, with the supremum over every positive common radius. Thus the argument is an operator theorem after the generic chain is accepted, not a kernel-only objective.

## 7. Full scalar comparison and a single exact candidate check

At x=(1-sigma)L,

    B0<=T+[sigma/(1-sigma)]n+1,
    A0+2sqrt(A0B0)=o(n)

for each fixed sigma>0. The complete new cost is at most T-c_raw*n+o(n); every fixed c<c_raw is therefore valid eventually. Neither the endpoint nor a common threshold across an adaptive continuous family is claimed.

`check_single_tuple.py` independently reconstructs **one** candidate directly from the generic formulas, with no import of the search formula and no grid traversal:

    a=73/8, s=293/16, d=1023/2^20,
    epsJ=epsD=epsS=epsT=1/1024, epsA=1/64.

It verifies every domain/sign condition, the rational budget, and c_raw>2^-16. It also records the frozen generic source hash and the candidate ledger hash, and checks a few rational constants used in this audit's relaxed-coordinate estimates. This is a single exact arithmetic corroboration, not an audit of the R052 enumeration, top-list completeness or resource claims; the numerical auditor owns those topics.

After this accepted generic chain and the single-tuple exact check, c=2^-16 lies strictly inside its admissible eventual coefficient interval. No finite-dimensional validity threshold is supplied. The R052 run remains a conditional numerical search in its original historical registration, which preceded the frozen generic proof; this audit does not rewrite that history.
