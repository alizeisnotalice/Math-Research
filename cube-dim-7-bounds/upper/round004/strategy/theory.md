# R004 Sol strategy: rare-shape diagnosis and one-dimensional pair-moment certification

2026-10-27. Independent clone of main at b2139fc1fe70fd55e7d409031256195c7d10e768. Only this `upper/round004/strategy/` directory is written. This agent executes **no numerical experiment, mathematical script, Monte Carlo, MATLAB, commit or push**. The statements below are analytical derivations and proposed finite acceptance problems, not computed outcomes. R001–R003 archived results and all other projects are unchanged.

## 1. Decision before allocating more CPU

The reported finite-shape minimum of H_(1/1024), about2.592e-6, is below the reported fixed-budget target g_req, about5.15e-6. This alone says neither that the true H minimum is that low nor that the full-G model fails. Two separate gaps remain:

1. A finite shape mean can miss the expectation; R002 already proved this for full G at/near S3.
2. Even the exact H_epsilon is only a **lower** approximation to G. An upper bound for H_epsilon at a point is not an upper bound for G there.

The reported block variation is a reason to diagnose the shape tails, not a quantitative upper bound on them. There should be no further full-G optimiser using the same inadequate cache.

The most useful additional allocation is a deterministic1D Gamma integral for an improved M bound (Sections4–6), together with a fixed raw-shape paired-cutoff diagnostic (Sections2–3). These answer different questions. An exact integral budget can improve the target even if the current apparent H minimum is real; expanding a point-only g experiment cannot repair a very loose I budget.

The proposed raw diagnostic must stop at its registered count. Saturating CPU means15 independent useful workers with BLAS/OpenMP threads1 and fixed scope; it does not justify extending an unresolved heavy-tail Monte Carlo run indefinitely.

## 2. Exact fixed point registry and discriminating observables

The parent has assigned the numerical agent15 processes:3 fixed points x5 independent seeds x2^20 new original Dirichlet shapes, with4 paired cutoffs. This strategy does not execute or enlarge that batch. Every process must evaluate its four cutoffs on **the same shape** and retain nonnegative adjacent-band differences. No point optimisation occurs.

### 2.1 Points and exact S3 status

The first point is the original R002 `gprobe-r002-start-013` terminal point

    x0=(1,1,1,1,0,1,1).

The read-only `starts/013.json` gives exactly these binary64 coordinates, and `diagnostics.json` first bundle names this source handle. With the R002 convention bit1=positive,

    x0=(v55+v63)/2,
    v55=(1,1,1,1,-1,1,1), v63=(1,1,1,1,1,1,1).

Thus x0 has minimal representation rank2, not rank1; its true G is infinite. It is a positive control for unresolved short-time/rare-shape contributions, not a low-G witness.

Let eta=1/1000003. For free-coordinate index j=1,...,6 define the fixed near point by

    y'_j=(1-7eta)y_j+j eta, x_near=(1,y'),
    y=(1,1,1,0,1,1).

Its free coordinates are1-(7-j)eta for j=1,2,3,5,6 and4eta for j=4. Together with the first coordinate1 these are seven distinct positive absolute values. Any combination of three body directions has at most four distinct absolute coordinate values, because the eight values +/-a+/-b+/-c form at most four opposite pairs. Therefore x_near is **exactly off S3**. This proves membership status, not a useful distance from S3.

The numerical agent supplied the R003 new Sobol point in exact normalized hex:

    ['0x1.0000000000000p+0',
     '0x1.cb8b8cf800000p-1',
     '0x1.ea76a78800000p-1',
     '0x1.d94dfdf800000p-1',
     '0x1.18baa38000000p-5',
     '0x1.7babca4800000p-1',
     '0x1.d943b29000000p-1'].

These normalized positive hex values are all distinct, so the same exact criterion proves offS3. Before the numerical batch the primary R003 source file, index and SHA256 must be locked by its owner; a chat quotation alone is not that provenance lock. An independent checker may convert each hex to its exact dyadic Fraction and test distinctness. No floating rank tolerance is needed.

An alternative rational generic point is(1,1/7,2/7,3/7,4/7,5/7,6/7), also offS3 by seven distinct absolute values. Calling it quantitatively far from S3 would require an exact/rational distance audit over all27216 rank3 spans; this strategy has not run one.

### 2.2 Required observables at the four cutoffs

Besides H estimates and block means, retain for each fixed point:

- paired band differences Phi_(epsilon_next)-Phi_epsilon, which are nonnegative;
- largest row share and a fixed number of highest-contribution row indices, without discarding any row;
- U=P55+P63 and a dyadic-bin decomposition of the contribution and adjacent-cutoff bands by U;
- counts and exact nominal Beta probability of bins, including empty bins;
- all nonfinite/SPD failures, with same-input higher precision rather than clipping/dropping.

Here U has **exact Beta(1,31)** law. This is more informative than recording only min(P_i): the two directions have a known geometric connection to the rank2 positive-control point. The near point has the same dominant geometric pair with a small residual. The Sobol point is not presumed to share this dominant pair; its U log is only a control, and a later selected triple must be fixed by exact geometry rather than an optimiser on the same data.

At x0 the exact H_epsilon diverges as epsilon decreases. A raw-cache plateau across cutoffs at that point is not model failure. At x_near/Sobol the exact G is finite, but convergence of H may be slow and point-dependent.

There is a useful qualitative fact: on compact cube Q, min_Q H_epsilon increases to min_Q G as epsilon decreases to0. To prove it, choose minimisers x_epsilon and a convergent subsequence. For every fixed delta>0 and epsilon<=delta, H_delta(x_epsilon)<=H_epsilon(x_epsilon). Continuity gives H_delta(x*)<=the limiting minima; increase delta down0 to get G(x*) below that limit. The converse is H_epsilon<=G. This proof does not give a numerical rate. True minimisers eventually avoid a fixed neighbourhood of S3, since G is infinite there and the compactness/continuous monotone-cover argument makes some H_epsilon exceed any fixed finite threshold on a neighbourhood of S3.

## 3. What bounded randomness can and cannot certify

R003 proved the uniform shape bound

    0<=X_epsilon(P,x)=Phi_epsilon(P,x)<=B_epsilon,
    B_epsilon=Gamma(71/2)/[6pi^(7/2)Gamma32·64^7epsilon^6].

For ideal independent samples, Hoeffding gives a simultaneous finite-point/cutoff confidence error h using N>=B_epsilon² log(2J/alpha)/(2h²). The dependence B_epsilon² proportional to epsilon^(-12) makes this raw range bound practically useless at the very small cutoffs. Apparent block stability or an empirical standard error does not replace a certified population variance bound. Unscrambled Sobol rows are not independent samples for this statement.

Clipping is useful for a **lower** diagnostic: for a prespecified cap b, X^b=min(X,b) is bounded by b and E X^b<=H_epsilon. A one-sided confidence lower bound for its mean can be much sharper than the raw B_epsilon bound. It remains a confidence statement under the sampling model, not an unconditional deterministic integral certificate. An upper bound for H requires paying E(X-b)_+; absence of sampled large rows cannot pay it.

If an independent analytic pth-moment upper bound K_p exists, p>1, then

    E(X-b)_+ <=[(p-1)^(p-1)/p^p] K_p b^(1-p).

This is the pointwise maximisation of(y-1)/y^p for y>=1. The trivial moment bound K2<=B_epsilon H_epsilon is usually too loose and must not be advertised as resolving the tail. Exact Beta-stratified bounds in Section7 give better bounds for selected low-U tails.

A deterministic positive-event contribution lower bound (Section8) can prove that a finite mean lies below a known part of the true expectation at a fixed point. That is an analytic diagnosis of missed contribution. It does not prove a global g. Conversely, a low H estimate cannot become an operator counterexample: a model exclusion requires a complete-G **upper** witness and a separately proved lower or exact budget expression (Section10).

## 4. Exact8-block decomposition and the priority M object

For each disjoint Hadamard block let

    R_g=sum_(i in g)Y_i~Gamma(4,1),
    p_g=(Y_i/R_g)_(i in g)~Dirichlet(1/2)^8.

These are independent within each block and across blocks. The R002 determinant identity yields

    det Sigma_g=16 R_g^(-7)/(product_i p_gi),
    T_g=(det Sigma_g)^(1/7)
       =16^(1/7)R_g^(-1)(product_i p_gi)^(-1/7).        (1)

Thus the parent's candidate radial expression is correct. Block determinant Minkowski gives

    M=E det(Sigma)^(-1/14)<=M_block:=E(sum_(g=1)^8 T_g)^(-1/2).          (2)

The T_g are iid strictly positive and nondegenerate. Direct8-variable AM–GM gave the R002 endpoint M_H. It need not be sharp, and the empirical value reported for actual M is not a certificate of how much slack comes from(2) versus AM–GM.

For any complex z with Re z>0, the exact Mellin moment is

    m(z)=E T^(-z)
        =2^(-4z/7) Gamma(4+z)/Gamma(4+8z/7)
                      ·[Gamma(1/2+z/7)/Gamma(1/2)]^8.                (3)

This follows from the positive Gamma/Dirichlet integral with a=z/7; absolute integrability is supplied by Re z>0. It is analytic in that half-plane and |m(a+it)|<=m(a). Also log m(a) is convex for real a>0 by Holder. These bounds are important for a rigorous quadrature tail; no asymptotic Gamma expansion is needed.

Block radial conditioning alone does **not** resolve the principal sparse rare event: R_g small scales all8 directions in a block together, whereas S3 singularity is produced by1–3 individual small Y_i, often inside otherwise ordinary blocks. One must retain the angular p_g boundary or use the selected-direction conditioning of Section7.

## 5. Two-block exact moment as a positive1D integral

Let U=T1+T2. For every real s>0 put a=s/2. The exact formula is

    q_s:=E U^(-s)
      =1/[2pi Gamma(s)] integral_R
            |Gamma(a+it)|² |m(a+it)|² dt
      =1/[pi Gamma(s)] integral_0^infinity f_s(t)dt,                  (4)
    f_s(t)=Gamma(a+it)Gamma(a-it)m(a+it)m(a-it)>=0 on R.

Here is a self-contained derivation. The beta integral under y=exp(r) gives

    Gamma(a+it)Gamma(a-it)/Gamma(2a)
      =integral_R exp(itr)[2cosh(r/2)]^(-2a)dr.

Fourier inversion therefore expresses[2cosh(r/2)]^(-s) as the Gamma-product inverse Fourier integral. Since

    (x+y)^(-s)=(xy)^(-s/2)
                    [2cosh((log x-log y)/2)]^(-s),

substitution and independence give(4). Absolute Fubini is legitimate: the expectation of x^(-a)y^(-a) is m(a)^2, and the Gamma-product integral is finite, for example by the recurrence bound below. Thus no formal Mellin contour manipulation is being assumed.

Now group the eight blocks into four independent pairs U_j and apply AM–GM only to those four numbers:

    M<=M_block<=M_pair:= (1/2)[q_(1/8)]^4.                            (5)

This is a **new deterministic bound expression** requiring only the positive1D integral(4), with s=1/8,a=1/16. It is analytically strictly smaller than the R002 expression:

    q_s < 2^(-s)m(s/2)^2,
    M_pair <8^(-1/2)m(1/16)^8=M_H.                                  (6)

The first inequality is the pair AM–GM inequality, strict almost surely because T1,T2 have a nondegenerate continuous distribution. Integrability was proved above. Strictness(6) alone does not quantify a useful improvement; this agent has not computed M_pair.

This still uses AM–GM at the four-pair stage and still pays determinant Minkowski in(2). It is not the exact M_block, nor the actual M. A recursive Mellin formula can exactly merge pairs into quadruples and then octuples, but nested complex integrals need new finite error registrations. The first1D object is the appropriate acceptance gate before paying that additional cost.

## 6. A fully specified finite certificate for the pair moment

The following is a proposed new finite batch for the main agent; **it is not executed here**.

Fix s=1/8,a=1/16, precision512 bits, T=32. For any positive integer N, Gamma recurrence and its defining integral give

    |Gamma(a+it)|
       <=Gamma(a+N)/product_(j=0)^(N-1)|a+j+it|
       <=Gamma(a+N)|t|^(-N).

Together with |m(a+it)|<=m(a), the exact positive tail obeys

    0<=q_s -1/[pi Gamma(s)] integral_0^T f_s(t)dt
       <=m(a)^2Gamma(a+N)^2/
                     [pi Gamma(s)(2N-1)T^(2N-1)].                   (7)

Choose N=16. This bound is elementary and valid for the entire uncomputed tail; no unproved numerical Gamma asymptotic is used.

For a deterministic finite-segment method, split[0,32] into512 equal panels, halfwidth h=1/32. Write the holomorphic extension

    f(z)=Gamma(a+iz)Gamma(a-iz)m(a+iz)m(a-iz).

It is analytic for |Im z|<a. Use a Cauchy radius r=3/64<a and degree128 Taylor polynomial on each panel about its real midpoint. For |Im z|<=r, the Gamma integral and moment absolute-value bounds give

    |f(z)|<=M_*:=Gamma(a-r)Gamma(a+r)m(a-r)m(a+r).                    (8)

To justify the endpoints in(8), Gamma(v)m(v) is log convex for v>0. Hence the symmetric product at a±y is maximised over |y|<=r at |y|=r. This supplies a single rigorous sup bound for **all** panels, independent of their real midpoint.

The unnormalised total Taylor remainder is bounded by

    E_quad<=32 M_* (h/r)^129/[1-h/r]
           =32 M_* (2/3)^129/[1-2/3].                               (9)

For Taylor coefficients b_j at a panel midpoint, integrate its polynomial exactly as

    sum_(k=0)^64 2h^(2k+1)b_(2k)/(2k+1).

The coefficients must be computed with outward complex-ball Gamma power series/automatic differentiation. Real and imaginary enclosures are retained; the mathematical coefficients are real on the real axis. An unchecked floating Gamma derivative routine is insufficient. If the installed runtime lacks the required certified series operation, report that implementation issue and use another verified method under the same error budget rather than claim a run.

Divide(9) by pi Gamma(s), add the positive tail(7), and assemble a rigorous q_s interval. The mathematical integrand is nonnegative, so intersect its lower bound with0 safely. Form M_pair=(1/2)q_s^4 and the complete original-G I endpoint

    I_pair=M_pair^(7/6)/(105·2^(25/6)pi^(13/12)).                      (10)

All real Gamma constants in(7)–(10) have exact rational inputs and are enclosed at512 bits. Check strict improvement against the exact R002 expression using outward endpoints; round the resulting upper endpoint upward to a declared rational grid. No g or C7 claim follows.

A four-panel pilot can verify the frozen series implementation and independent recurrence/positive-integral identities. The full frozen batch then assigns512 panels round-robin over15 processes (34 or35 panels each), BLAS/OpenMP1, with immutable code/config hashes. There is no reason for unbounded MC or a63D partition here. Stop after512 panels/tail/assembly and main review; do not automatically start quadruple merging.

## 7. Low-dimensional conditioning aimed at the sparse rare carrier

Fix k<=3 independent dangerous directions T and let R denote the other64-k directions. The Dirichlet law splits exactly as

    P_i=U q_i for i in T,
    P_i=(1-U)r_i for i in R,
    U~Beta(k/2,(64-k)/2),
    q~Dirichlet(1/2)^k, r~Dirichlet(1/2)^(64-k),

with U,q,r mutually independent. Put C=A_R(r), V_T for selected columns and D_q=diag(q_i). Then

    A(U)=C/(1-U)+V_T D_q^(-1)V_T^T/U.

For fixed r,q, evaluating the U integral requires only kxk Woodbury data. Let G_T=V_T^TC^(-1)V_T, b=V_T^TC^(-1)x, q0=x^TC^(-1)x, tau=U/(1-U). Then

    det A=(1-U)^(-7)det C·det(I_k+tau^(-1)D_q^(-1)G_T),
    Q_A=(1-U)[q0-b^T(G_T+tau D_q)^(-1)b].                           (11)

A stable evaluation should separate d0=q0-b^TG_T^(-1)b>=0 from the nonnegative rank-update term, rather than subtract nearly equal numbers at U near0. If x=V_T alpha, the leading expression is Q_A~U alpha^TD_q alpha as U decreases to0. One can also compute that case directly with a kxk positive solve.

For x in spanT and fixed interior r,q,

    det(A)^(-1/2) asymptotic const·U^(k/2),
    Q_A asymptotic const·U,

so the full-G shape integrand times the U density behaves as U^(k-4): divergent for k<=3. The cutoff changes the low-U product to O(U^(k-1)epsilon^(-6)), which is integrable. Its transition occurs at U on the scale epsilon² divided by the positive geometric coefficient. This is precisely why a standard quadrature with no near-zero stratification, or a small raw sample, can miss the relevant contribution.

For the rank2 control T={55,63}, U~Beta(1,31) and q1~Beta(1/2,1/2). The transformation q1=sin²theta makes theta uniform on[0,pi/2]. Hence integrating U and theta is only a2D deterministic conditional integral. A1D U conditional integral already removes its most dangerous radial boundary; q can remain as an outer variable. For k=3, q1=v², q2=(1-v²)cos²theta, q3=(1-v²)sin²theta with independent v uniform[0,1], theta uniform[0,pi/2]. Thus at most3 well-specified variables replace a microscopic63D box. The high-dimensional remainder r remains an expectation; conditioning does not secretly certify that outer expectation.

A finite diagnostic pilot can reuse a **fixed prefix** of the new independent r,q rows (no new draw), integrate U on a declared dyadic partition down to a fixed endpoint, and enclose the omitted U tail by the bounds below. Low-dimensional deterministic conditional integration reduces variance by the law of total variance; it must not be called a deterministic full expectation certificate when r is still sampled.

A useful tail bound is available without any tiny box for r. Since C>=a0 I, a0=64-7k, the determinant's selected-rank term gives

    det A >=a0^(7-k)(1-U)^(-(7-k))U^(-k)
                         det(V_T^TV_T)/product_i q_i.

Together with the R003 profile amplitude, average first over q and r to obtain

    E[Phi_epsilon 1_(U<=eta)]
      <= K0/[6epsilon^6 a0^((7-k)/2)sqrt(det(V_T^TV_T))]
        ·Gamma(k/2)/[Gamma(k)pi^(k/2)]
        ·B_eta(k,(71-2k)/2)/B(k/2,(64-k)/2).                         (12)

Here B_eta is the unregularized incomplete Beta. The q factor is E sqrt(product q_i). Formula(12) decays as eta^k, stronger than a crude B_epsilon times the Beta tail of order eta^(k/2). Every constant is explicit. Analogous second-moment bounds follow by squaring the determinant envelope, with q factor Gamma(k/2)/[2^kGamma(3k/2)] and incomplete-Beta parameters3k/2 and(78-3k)/2. These are strict population bounds for the selected tail, not sample-variance estimates.

A generic eight-Hadamard-block decomposition instead writes group masses beta~Dirichlet(4,...,4), independent of eight Dirichlet(1/2)^8 angles. Conditional group-mass integration is7D and its stick-breaking Beta(4,4(8-j)) box masses are polynomial/rational at rational endpoints. It is a useful secondary scheme, but it is larger than the selected-k conditioning and does not remove the single-coordinate angular boundary by itself.

## 8. Deterministic positive-mass contribution without a63D microbox

Here is a separate exact **point lower-bound** construction using the raw independent S_i. It can make a missed contribution analytic, although its first crude constants may be weak.

Choose independent selected columns T, |T|=k<=3, and let N=64-k. Let Sigma0=sum_(i not in T) S_i v_iv_i^T. For B,L>0 consider the remainder event

    E={sum_(i not in T) S_i<=L and max_(i not in T)Y_i<=B}.

The sum has the exact stable law N²S, so

    Pr(E)>=p_E:=max(0, erfc(N/(2sqrt(L)))-N erfc(sqrt(B))).          (13)

This is a lower bound from a union inequality; the two remainder constraints are not falsely assumed independent. On E,

    Sigma0>=aI, a=(64-7k)/(4B), trace Sigma0<=7L.                    (14)

The first follows from the remainder Gram64I-sum_T vv^T>=(64-7k)I and S_i>=1/(4B). E is independent of the selected S_i.

For D=diag(selected S_i), the determinant has a positive multiaffine expansion in those k variables. For J subsetT its coefficient is det(V_J^TV_J) times the determinant of the compression of Sigma0 to span(V_J)^perp. By AM-GM on that compression and trace<=7L,

    det(Sigma0+V_T D V_T^T)<=D_+(S),
    D_+(S)=sum_(J subsetT) det(V_J^TV_J)
             [7L/(7-|J|)]^(7-|J|) product_(j in J)S_j.              (15)

For J empty the term is L^7. There are at most8 terms, and all coefficients are explicit Gram determinants. This is substantially sharper than replacing Sigma0 by7L I and paying its determinant in every dimension.

For x=V_T alpha+w with w orthogonal to spanT, Loewner monotonicity gives

    x^T(Sigma0+V_TDV_T^T)^(-1)x
       <=Q_+(S)=|w|²/a+alpha^T(D+a(V_T^TV_T)^(-1))^(-1)alpha.       (16)

For fixed covariance Sigma the truncated Gaussian-time integral is

    C det(Sigma)^(-1/2)Q^(-3)P_gamma(3,Q/(4epsilon²)),
    C=1/(2pi^(7/2)),

with its continuous Q=0 extension. Its scalar Q profile is decreasing. Thus on E, replacing det by D_+ and Q by Q_+ gives a valid lower integrand. Integrating selected k Gamma(1/2) variables and multiplying p_E gives a deterministic k-dimensional lower bound for H_epsilon(x). At k<=3 this is a concrete interval quadrature, not a63-dimensional tiny event.

An even simpler finite lower certificate uses selected S_i in[s_lo/epsilon²,s_hi/epsilon²], with0<s_lo<s_hi. If w=0,

    Q_+<=epsilon² sum_i alpha_i²/s_lo=:epsilon² q_bar.

Using only t in[epsilon,2epsilon] yields

    H_epsilon(x)>=p_E
      ·[Pr{Y in[epsilon²/(4s_hi),epsilon²/(4s_lo)]}]^k
      ·(4pi)^(-7/2)D_+(s_hi/epsilon²,...,s_hi/epsilon²)^(-1/2)
      ·exp(-q_bar/4)·(1-2^(-6))/(6epsilon^6).                       (17)

All probabilities are Gamma1/2 CDF differences, i.e. erf differences, and every term is explicit positive. Proposed finite positive-control choices are B=8,L=N²,s_lo=1,s_hi=2 and the four already registered cutoffs; they have not been evaluated here. Exact Gram determinants for the selected pair make(17) a small scalar Arb calculation. Failure to produce a useful lower contribution rejects this particular coarse event, not the model.

For k=1,2 the expression scales as epsilon^(2k-6) as epsilon goes to0; for k=3 a single selected scale band gives a constant, and summing disjoint dyadic S bands up to the cutoff produces the logarithmic growth. A near-point w!=0 has the extra term |w|²/(4a epsilon²), explaining why the positive-control and nearby finite-potential point separate only at a sufficiently small cutoff. A deterministic lower contribution exceeding a reported empirical point value proves that value is missing known expectation mass; it still supplies no global face minimum.

## 9. Separate the four sources of M-budget slack on the same P rows

No extra random shape is needed to diagnose where the empirical M versus M_H gap comes from. For global P let

    B_g(P)=sum_(i in block g) v_iv_i^T/P_i,
    r_g=det B_g(P)^(1/7), A=sum_g B_g,
    D0=2Gamma(65/2)/Gamma32.

The following per-row quantities satisfy exact inequalities:

    Z0=det(A)^(-1/14)
       <=Z1=(sum_g r_g)^(-1/2)
       <=Z2=(1/2) product_(j=1)^4(r_(2j-1)+r_(2j))^(-1/8)
       <=Z3=8^(-1/2)product_(g=1)^8 r_g^(-1/16).                   (18)

D0 E Z0 is actual M; D0 E Z1 is exact determinant-Minkowski block budget; D0 E Z2=M_pair; D0 E Z3=M_H. Thus means of their paired differences isolate actual-M to determinant-Minkowski slack, pair-stage AM–GM slack, and the remaining pairwise AM–GM slack. These are diagnostics unless independently certified, but far more informative than comparing only Mhat and M_H.

All Z_i<=1/64. To see this for Z3, write group P mass beta_g and internal p_g. The tight8-column frame AM–GM determinant bound gives r_g>=64/beta_g. Then product beta_g<=8^(-8), so Z3<=1/64; the others follow(18). Hence these M diagnostics have the modest uniform bound D0/64, unlike the epsilon^(-6) H range. Bounded-confidence estimates can genuinely distinguish large budget gaps at finite N, while still not being unconditional exact expectations.

If the15-process code was already frozen, do not silently add eight7x7 determinants per row. Register a fixed prefix of already new rows for(18) as a separate finite diagnostic. Store exact partition identity, Gram/tight checks and original row indices; no new MC is needed. Never conclude that M_block≈M merely because both are upper bounds.

## 10. Exact model/budget exclusion gates, with correct quantifiers

A strict model exclusion needs a **complete-G upper witness** at a fixed offS3 point x*. A certified H upper witness alone is insufficient. One possible route is H_epsilon(x*) upper plus a certified short-time upper remainder; R001 supplies a finite O(epsilon²) local bound offS3, although its crude constant may be unusable. A sharper conditional integration must pay every unintegrated positive tail, not merely report a low raw mean.

There are three distinct exclusion claims:

1. **Specified analytic budget only.** Let J be the exact deterministic integral-budget expression, e.g. R002 I_H or the new I_pair, not the unknown true I. If g<=G(x*)<=g_U and a certified lower endpoint J_L of that deterministic expression satisfies

       7^7 J_L^6 >=128^6 g_U^7,

   then even the best true g cannot make the bound using that specified J beat128. J_L is a lower bound on the **budget expression**, never a lower bound on actual I. This rejects that budget, not the frame and not C7.

2. **The full-boundary Minkowski family with exact actual M.** Put F=1/(105·2^(25/6)pi^(13/12)) and kappa=7(F/128)^(6/7). Its best exact-M budget fails if a rigorous M_L<=actual M and g_U<=kappa M_L. A general rigorous M_L is1/(32sqrt(pi)): Hadamard gives det Sigma<=(sum S)^7 and sum S has law64²S. A tighter certified M_L would improve this gate. For determinant-block refinements only, Jensen also gives M_block>=(8E T)^(-1/2), where E T=m(-1) is finite and its Gamma formula extends to Re z>-7/2. Excluding the block-refined family by that lower budget does not exclude the original exact-M route.

3. **The original positive-part frame model itself.** This requires actual I_L<=I as well as g_U, and the same exact power comparison7^7 I_L^6>=128^6 g_U^7. A certified integral of H_epsilon^(7/6) is a valid possible I_L. A numerical estimate of it is not. Even a successful test excludes this particular equal-weight/common-time positive-part model as a way to beat128; it is not a counterexample to the operator weak bound or a universal impossibility theorem for other kernels/frames.

No such upper-G witness or model exclusion has been certified in this strategy task. The empirical M report motivates a tighter I budget; it does not authorise pretending the loose M_H endpoint is actual M.

## 11. Recommended next acceptance order

1. Finish the already frozen15-process point/paired-cutoff diagnostic, exact3-point status and source locks. No optimiser or adaptive sample enlargement.
2. Prioritise the1D pair-moment certificate(4)–(10): four-panel implementation pilot, then the frozen512 panels+explicit tail. This is the direct zero-MC path to a tighter complete I bound and uses15 CPUs on independent verifiable panels.
3. In parallel, evaluate the finite positive-event expression(17) for the rank2 control at the four fixed cutoffs. Its certificate can expose a known missing contribution. If the coarse event is too weak, report its explicit gap; do not grow a63D box scheme.
4. Only then consider fixed-prefix selected-U conditional quadrature or the slack observables(18), under separate finite registration. Exact conditional integrals can diagnose the rare carrier with many fewer new outer rows. Do not call an outer MC expectation deterministic.
5. Use Section10's exact threshold tests to decide whether the remaining g problem is genuinely competitive. If a budget is excluded, identify precisely which budget/model layer is excluded. If unresolved, keep that status rather than infer failure from the low H means.

All next-batch counts, precisions, tails and stop conditions above are proposals for the main owner; this agent has not launched them. The new analytical pair bound, low-dimensional rare-carrier factorisation and deterministic positive-event formulas are ready for independent review. No numerical result or new g/C7 bound is claimed.
