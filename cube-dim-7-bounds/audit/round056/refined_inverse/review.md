# R056 refined inverse mechanism: independent full-chain review

2026-10-28. **ACCEPT at the corrected final source.** For every integer n with natural L=log n>=2^20, the proof establishes

    x=93L/100, t=n/x,
    C_n<=F(A_n(t),B_n(t)-gamma n)
       <n log n-n log log n-n/40,

for the exact gamma below. This is an upper bound for the original Lebesgue centred axis-parallel cube maximal operator, with general L1 inputs and all positive common radii. Leading coefficient 1 remains. No production diagnostic, parameter optimum, numerical-list identity or C_n lower bound is accepted by this conclusion. The main coordinator retains final integration/proof responsibility.

Final audited source: `/Users/zhengzhihao/Desktop/T/cube-r030-geometry/cube-dim-7-bounds/upper/round056/refined_inverse/PROOF.md`, SHA-256 `6986ed9ec48b46f3136a9a56722d0092db4f2942d52e35c9f4de37eb279f399a`.

## 1. Sources, division of review and corrected identity

The initial full source `5d639f31cb7c62d05a09612d63fe22fc1af408cc5f311eed72d803c7e0c12eb3` was read in full. This auditor found one genuine local probability-wording error: the factor one half for the original target cap is correct after averaging positive/negative maximum faces, but not on each signed face as claimed immediately afterwards. On a positive pinned face the cap fraction is ell_t(r)^(n-1); on a negative one it is zero. The author corrected the source to the stronger uniform per-face bound (3/4)^(n-1), and removed the half from its guard. The bound (3/4)^74<=epsA/4 still has large strict margin, so no theorem budget or threshold changed.

That correction yielded intermediate source `623b7b15967e1437e45c641d171c84a50b5847558ccade39d7913f6be78028b8`. The final source additionally makes the full source band, bounded b extension, cap-restricted kernel/degrees and actual source-radius relation explicit. Both deltas were inspected. An independent line diff against the original confirms only those per-face statements, the new definitions, and the source-item-radius explanation changed. This auditor modified no author file.

Assigned local inputs were also read and matched to the actual construction:

* roots `audit/round056/one_coordinate/LEMMA.md`, final SHA `97430af9980c2de5d105bc4ea86667573c5c28e0ff40a92f1c613cd95a243eec`: exact J, expanded local domain, baseline moments, a.e./AC bridge and exact free-source phi moments;
* soft `audit/round056/refined_kernel/mechanism_review.md`, SHA `f5615fc233a3a624684ac5c06df882b1f0c7395c7fa5b3660a4aa9a630de7283`: endpoint face calibration, column factor, per-face inverse Chebyshev and post-cap payment;
* parameters `upper/round056/parameter_gates/PROOF.md`, actual SHA `e0aa014154ee7ce8e4ce558a9c48cf76fbc93d4cf603c826f927901086e4f077`: separate exact budget/gate proposal, not full operator acceptance by itself.

The assigned reviews alone do not certify the operator theorem. This review checks every interface, actual-history/finite-a/L1 passage, explicit all-tail start and full scalar cost, retaining its own judgment. No local-J baseline program, author 29/49-gate checker, old parameter grid, production batch, numerical exponentials, n=7 work, shared MATLAB or Git operation was run. Only this independent clone's `audit/round056/refined_inverse/` directory was written.

## 2. Geometry at d=3/10 and the correct baseline

The proof does not apply the R051 d<=1/16 theorem outside its domain. It uses a separately derived local lemma valid for d<=1/3, R>=1, g<=R/2 and rho>=4r. With h=2^-24/(Ld), source q in [31/40,37/40] and trimmed target radii

    rlo=43/40+2^-20, rhi=49/40-2^-20,

the matching interval is [r-d-h,r-d+h]. The **baseline is its upper endpoint qb=r-d+h**, where R=q+d-h=r. The centre r-d would not give the exact E J=1 cancellation and is not used for that purpose.

At L0=2^20, N=L^8/8!<=n is already greater than 75. Increasing lower bounds N/L and qmin*N/L give t>=4 and rho>=4rhi. The guards h<=2^-20 and 1/(xN^2)<=h/2 imply epsilon=zeta*t<=h/2, beta<2, kappa>1, 0<g=d-h+epsilon*q<=d, h<d/8 and g<R/2. The full matching interval fits inside the source band: its extreme endpoints are qmin+2^-20-h and qmax-2^-20+h. Its ordinary width is exactly 2h even at trimmed target endpoints.

The reference radius is at least rlo-2h>10749/10000>1, on every target matching interval as well as the full source band. All local values have q>=1/2, rho>=4r, a0=(1-epsilon)q>0 and T=(|w|-a0)_+<=g+2h<=5d/4. Thus the expanded local lemma's hypotheses are actually met. These margins strengthen for larger L; no residual eventual N is being silently used at the fixed endpoint.

## 3. Local J, pins, AC and target probability

The exact formula is J=1+[F(R,g)+F(R,(g-m)_+)-m]/(2rho), m=(|w|-a0)_+, with F(R,s)=-2R log(1-s/(2R))-s. It covers deletion beyond the tail into the central region. The local proof establishes F_s<=d, |F_R|<=d^2 and F<=3d^2/10 in the enlarged domain. Its capture-onset kink is handled by local Lipschitz continuity and absolute continuity, so the derivative estimate is a.e. and may be integrated to give a statement for **all** matched q. No old everywhere-differentiable assertion or small-d estimate is transplanted.

At R=r, a free W uniform on [-r,r] satisfies E J=1 exactly by Fubini. The -6/12 log mean/variance constants follow from the elementary second moment and log inequalities at R>=1,g<=R/2. The pinned target coordinate is deterministic, with log in [-g/rho,0]; it contributes a mean term and no variance.

For the target log sum the uniform mean/variance upper bounds ML,VL in the source are correct: x<=L and rho>=n*qmin/L give the factors 6d^3 L^2/(qmin^2 N), dL/(qmin N), and 12d^3 L^2/(qmin^2 N). The stated base conditions ML<=1/512 and VL<=2^-26 imply the log-event failure <=1/256=epsA/4 at tolerance 1/256.

For each endpoint q_-,q_+, the target pin has N_i=1, since a0<r. The other n-1 indicators are independent Bernoulli of parameter a0/r. Hence E T_N=n+(n-1)a0/r and Var T_N<=n/4; threshold e*n/4 gives each endpoint bad probability <=4/(e^2 n). The original corner bad fraction is now bounded correctly **on each face** by (3/4)^(n-1). Since t>=4 and rhi<49/40, ell_t(r)<3/4, and n>=75 supplies the corrected stronger guard. No independence among the four good events is assumed.

Our independent compressed check evaluates the direct log Chebyshev bound plus both endpoint probabilities plus the stronger per-face cap bound and finds their sum <epsA at the base. Uniform decreasing bounds extend this to every face/radius on the entire tail. Thus lambda(A)>=(1-epsA)n log(rhi/rlo)>=mA*n. The four-term alternating lower polynomial for log(1+zA), 0<zA<1, is a valid lower bound and positive. It is not confused with the exact logarithmic mass.

For d=.3, K_w<=n<=4nd is automatic. Integrating the local derivative over length 2h yields 1280Ldh=5/65536. Baseline 1/256 plus this change gives the precise P0=65275/65536 product enclosure. Individual J bounds are supplied separately by |J-1|<=2d/rho; they are uniform over every allowed coordinate and radius, not an n-coordinate failure union or a product of n coarse bounds.

## 4. Endpoint count calibration, exact rows and the improved column

The source face count and pin factor are calibrated jointly. T_N(q) and ell_*(q) increase with q. Applying the lower count at q_- and upper denominator at q_+ gives loss at most 2h+1/n+e/4. The 1/n term correctly retains the n-1 free-coordinate correction. For the reverse pairing,

    a_+/r-a_-/(r-2h)
       =2h(1-epsilon)(d-h)/[r(r-2h)]<=2hd.

Thus the upper loss is at most 2hd+e/4. The common base guard 3h+1/n<=e/2 controls both for every intermediate q. These are two events on the same fixed target; no unrelated worst endpoint values are summed.

The complete inverse-source angular coefficient is

    U=(product J)/(2ell_*) sum_i N_i/J_i.

Its pin factor is retained. Combining the calibrated T_N/(2n ell_*) with individual J and product bounds gives exactly Ulo=P0(1-e)/(1+e) and Uhi=(1+e)/[P0(1-e)]. There is no lost factor one half in Ulo and no unjustified independence between target count and product events.

The new column factor 2 is pointwise. For fixed source z, a_i=kappa*rho-|z_i|>0 and target normalized overlap is ell_i=min(1,(1+a_i/r)/2). If a_i>=r, ell_i=1 with at most two target faces; if 0<a_i<r, exactly one face is captured and ell_i>1/2. At equality the strict capture convention can only reduce the count. Hence N_target,i/(2ell_i)<=1, and the exact shell coefficient is <=n product ell_i. Each ell_i decreases as r increases from the lower reference R. Matching target log width is <=2h/R, proving C_raw<=2nhG/R. The bounded soft factor gives column<=2k. This is an ambient maximum, not an average face-count estimate.

## 5. Per-source-face inverse Chebyshev before cap restriction

The two full-source free-coordinate identities are exact:

    E_full[phi/ell]=R Qpin^2/rho,
    E_full[phi^2/ell]=2R Qpin^3/(3rho),
    phi=-log ell, Qpin=-log ell_*.

Double-tail normalization has already been included. On one fixed captured source face, free inverse-capture laws have density capture/(ell J_i) with J_i>=1-e. Removing the nonnegative capture indicator gives a legitimate upper bound for each moment. Free coordinates are independent **on that fixed full source face**; the source pin contributes Qpin deterministically. Thus its mean is <=Qpin+x R Qpin^2/[q(1-e)] and variance <=x*2R Qpin^3/[3q(1-e)].

R/q<=43/31, g<=.3 and R>=10749/10000 give Qpin<=log(3583/3083). The exact cubic upper polynomial Q for log(1+500/3083) has Q<1/6; the mean and variance rate bounds <4/125 and <1/300 are valid. At x>=1000, Q/x<=1/2000; chi-4/125-1/2000=1/400. Chebyshev on **each face separately** therefore gives tail <=1600/(3x)<=e. The mixture uses these uniform probability bounds, not a false equality between total variance and averaged within-face variance.

The full-source inverse law is used before either cap conditioning or the bounded soft weight changes the distribution. On G>=theta, soft mass is theta times inverse mass, and elsewhere it is nonnegative. This retains at least (1-e)Ulo*n*theta on every matched shell; the upper is Uhi*n*theta. Full matching width is exactly 2h, so raw weighted rows lie between (1-e)2Ulo*k/qmax and 2Uhi*k/qmin.

The final source explicitly defines b on the **whole** source band, including outside D, and only extends it as zero off that band. The final kernel then restricts the source to D. This correctly preserves all full-source laws used above. A and P are disjoint, all kernel indicators/weights are Borel and 0<=K<=1, and the displayed degrees are the integrals of that same kernel.

## 6. Refined cap, tau, degree and explicit infinite tails

The source cap needs the new exponent. At n>=75 and q<=37/40,

    (1-1/n)(1-q)/(2q)
       >=(74/75)*(3/74)=1/25.

The shell cap coefficient is <=(n/2)exp(-x/25), and the source log-band width is <=d. Hence its unweighted full mass is <=nd exp(-x/25). Since b<=1 this bounds removal from every weighted row. Dividing by k pays the precise gap 1/25-chi=1/200. The old Ld/32 exponent would fail here and is not used.

Global Taylor order eight yields the source's cap bound proportional to L^-7. For tau, A<=n exp(-sigma L), B<=nL and theta=exp(-chi x) give the exact exponent sigma/2-chi(1-sigma)=49/20000. Taylor order eight yields a bound proportional to L^(-13/2). Our new exact checks verify both at L0 below e. These full-domain power bounds cover every larger L without numerical exponential evaluation or a finite grid.

The final row lower subtracts cap loss e*k; its excess after actual threshold payment subtracts a separate e*k. That gives b_*=(1-e)2Ulo/qmax-2e>0. The upper row remains 2Uhi*k/qmin<21k/8, while column<=2k. Thus the common positive degree denominator M*k=(21/8)k is valid. Target mass yields missing measure >=mA*b_*n/M; keeping gamma=(1-e)mA*b_*/M leaves a strict epsD margin.

Other infinite tails are fully quantitative: N and N/L increase; baseline means/variances are sums of L^-6,L^-7; count probabilities are O(L^-8); individual J error is O(L^-7); h and count-calibration drifts decrease; the epsilon/h bound is O(L^-16); per-face retention bounds decrease as x^-1. The target cap only improves as n increases. Lower geometric margins improve with decreasing h. There is no residual O_d term or unquantified eventual start at L0.

## 7. Same actual history, finite dilation and general L1/all radii

B>8n follows directly from the original cap integral bound for x in [L/2,L] at this L; B<=nL and gamma<1 make B-gamma*n positive. No small-d family theorem is used for that scalar estimate. The final actual measure consists of the original uncaptured positive-contribution predecessor increments; it obeys mu_a<=lambda_a. On both new supports lambda_a=alpha^-n lambda exactly, by D subset H_a.

For an actual target item j the uncaptured-centre inequality gives ccap*r_j<alpha+r<3, with ccap>1. For source increment points z=y/(mathfrak a+1), the original common support and increment inclusion imply ||z||<r_i; thus r_i>=m(z)>=rho>4. Every source precedes the **same** target j. Since x>1, t<n and every mathfrak a>=8n^4 has mathfrak a>2t. With beta<2,

    alpha beta-kappa-2/(mathfrak a+1)
       =n^-3-(beta+2)/(mathfrak a+1)
       >n^-3-4/(8n^4)>0.

Increment-to-centre displacement then gives true capture by the original dilated source cube. Bounded b<=1 pays the actual weighted row by that target's original eligibility tau; no inverse-theta factor is charged and no distribution/concentration assumption on mu is made.

For c_a=alpha^-n>=1, ambient energy scales by c_a^2 while degrees and target mass scale by c_a. The energy identity uses a single missing measure and disjoint supports, so its finite-a bound is at least [c_a E-tau lambda(A)]/(M*k), which is >=the unscaled bound. Missing cap mass is therefore >=mA*b_*n/M>gamma*n for every such original finite history/dilation. This gives the stated uniform finite cover/load inequality, including captured tau and terminal self-weight mathfrak a^-n.

Holding n,t,tau fixed, only these scalar inequalities are passed to mathfrak a->infinity. The original compact strict-level finite cover uses the same general L1 input and the original mathfrak-a-dilated witness cubes; inner regularity returns the full level set. No selector/history convergence, new input class, compact-support requirement or bounded-scale operator is substituted. The first operator inequality is thus proved on the claimed full domain.

## 8. Full scalar cost and exact budget

The rational lower target mass and chosen deficit reconstruct as

    mA=259465051078193863093554363/
       2018111068113594386771608802,
    gamma=172531397335235056462877787662549380239/
          1643616601687301779182642153860479582208.

The exact strict slack gamma-7/93-1/40 equals

    3593475707213855800847480839340601290687/
    764281719784595327319928601545123005726720 >0.

This mA is a valid polynomial lower bound, not the logarithmic mass itself. The cost includes the full buffer -log(1-sigma)<=sigma/(1-sigma)=7/93. Also sigma/2=chi, so the positive root-cross term is <=2n sqrt(L) exp(-chi L). Combining A, that cross term and 1/n gives exactly the stated Taylor-eight P(L), with powers L^-8,L^(-15/2),L^-8. Its exact base value is below the slack above and decreases on the entire tail. Thus the final full-root inequality is strictly below T-n/40, rather than being inferred just from gamma or from a cap-component reduction.

## 9. New independent execution and final lock

`microcheck.py` independently reconstructs the fixed budget and runs 20 compressed exact certificates using Python/Fraction. It does not import any author helper. Its target check uses an aggregate four-event union, and its inverse-tail check also uses direct exact mean/variance rates; these differ from the registered author gate layout. It additionally checks trimmed/local geometry, individual J, both endpoint calibrations, the pointwise degree budget, exact x/25 and tau exponents, both Taylor-eight margins, finite-dilation positivity and the complete scalar slack. Final execution PASS, PID 90491. All exact values and actual source hashes are in `microcheck.json`.

One earlier invocation reached the source-identity assertions after all mathematical checks but stopped because the roots lemma had concurrently acquired its final source-binding section. That append was read; it changes no local estimate. The hash lock was updated and the same script executed successfully. Neither invocation ran the author's guards or local baseline verifier. The source-change history is recorded in Section 1 so that initial and final identities are not conflated.

Final full-source SHA remains `6986ed9ec48b46f3136a9a56722d0092db4f2942d52e35c9f4de37eb279f399a`. The corrected per-face cap statement, full-band-before-cap kernel definition, and actual-prefix explanation are included in this acceptance. No outstanding mathematical correction remains. Terra kernel diagnostics and proposed later-round candidates are not used as proof or included in this review's scope.
