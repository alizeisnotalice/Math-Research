# R016: an origin-density seventh-moment bound for the complete I budget

2026-10-27. Sol theory only. No numerical scalar/kernel/phase/MC/MATLAB run by this owner; the primary scalar pilot is separate source evidence, not this owner's calculation. Fixed original7D body64, unit weights and common Cauchy time. Lebesgue measure, centered axis-parallel cubes, all common radii and generalL1 retain their accepted original comparison interface.

## 1. The original Gaussian mixture and a pointwise Tonelli identity

Let S_i>0 be independent positive stable variables with E exp(-lambda S_i)=exp(-sqrt(lambda)), and Sigma=sum_i S_i v_i v_i^T. The64columns span R7, so Sigma is SPD almost surely. Averaging exp(-t² xi^T Sigma xi) gives exp(-t psi(xi)); thus this is the original common-time Cauchy process. The conditional density at time t is

    (4pi t²)^(-7/2)det(Sigma)^(-1/2)
          exp[-x^T Sigma^(-1)x/(4t²)].

In particular the origin identity, including its precise pointwise version, follows directly by positive Tonelli in frequency:

    p1(0)=(2pi)^(-7) integral_R7 exp(-psi(xi))dxi
         =E[(2pi)^(-7) integral_R7 exp(-xi^T Sigma xi)dxi]
         =(4pi)^(-7/2)E det(Sigma)^(-1/2).                    (1)

This avoids relying on an a.e. equality of two density versions at a selected point. The integrands in the second line are nonnegative real. A Euclidean floor psi>=c|xi|, c>0, proves p1(0)<infinity by radial integration, so the expectation in(1) is genuinely finite.

Write Z=det(Sigma)^(-1/14), M=EZ. Then Z^7=det(Sigma)^(-1/2). Hölder against the constant1 on the probability space gives

    M<= (EZ^7)^(1/7)=2sqrt(pi) p1(0)^(1/7).                  (2)

Equivalently apply concavity of x^(1/7) to the integrable positive variable Z^7. The direction is an upper bound. No expectation of Sigma is used or required. p1(0) is a finite time-one density; it is not G(0). Indeed p_t(0)=t^-7p1(0), so the full potential G(0)=integral_0^infinity p_t(0)dt is infinite.

## 2. Insert only the new M bound into the accepted whole-boundary interface

For q=7/6, the accepted R001 Gaussian argument on the entire six-dimensional cube boundary yields

    I=integral_[0,1]^6 G(1,y)^q dy
       <=M^(7/6)/[105*2^(25/6)*pi^(13/12)].                   (3)

For clarity, its constants arise from the actual all-time potential factor1/(2pi^(7/2)), the cubical ray identity integral_boundary (x^T Sigma^-1x)^(-7/2)dσ=(16pi³/15)sqrt(detSigma), whole-boundary Lq Minkowski, and only then the original-process symmetry factor7*2^7=896. No individual conditional Gaussian is assumed to have cube coordinate symmetry. Positive truncations and monotone convergence justify Minkowski without presupposing the final I is finite.

Substituting(2), the powers simplify exactly:

    (2sqrt(pi))^(7/6)/[105*2^(25/6)*pi^(13/12)]
      =1/(840sqrt(pi)).

Consequently the new complete budget is

    I<=B_origin:=p1(0)^(1/6)/(840sqrt(pi)).                   (4)

This is a full original-G, complete-time, whole-face budget. It does not use a time cutoff, a pointwise minimum, a shape sample or a new operator comparison. Only the numerator estimate has changed.

## 3. Exact fan volume gives p1(0)

For each accepted ordered simplex S let alpha_S=|detW_S|/product q_i, and let A_sum=sum_212 alpha_S. All alpha values are positive exact rationals. On its cone, direct integration of exp(-psi) gives alpha_S. The645120signed-permutation ordered cones cover frequency space with disjoint interiors. Hence

    p1(0)=645120 A_sum/(2pi)^7=5040 A_sum/pi^7.               (5)

The scale is primitive frequency xi=Mw, q=psi(xi), detM=1. Neither the u-simplex volume determinant nor its normalization denominators should be inserted a second time. The accepted R011 exact sum is

    A_sum=907157967345860903220972676779187
           /28127466839029475035590994122592223232000000000.

Combining(4)-(5),

    I<=B_fan=(5040 A_sum)^(1/6)/(840pi^(5/3)).                (6)

The formula is proved conditional only on the accepted exact fan certificate, not on any sampled generic kernel values. It requires no new212simplex or full-phase run. The primary thread's saved512-bit scalar pilot reports approximately p1(0)=5.381876473321905e-14, M bound0.04508416552889 and B_fan=4.12700249882943e-6. These are provenance-labelled pilot displays; this owner did not execute them and they are not substituted for a separate strict scalar certificate.

## 4. A sharp-floor calibration independent of A_sum

The accepted spectral-unit-ball certificate and its new R011 scalar functional prove psi(xi)>=c|xi| with c=32sqrt2. This floor does not require the R010 exact fan volume. Since area(S6)=16pi³/15 and Gamma7=720,

    p1(0)<=(2pi)^(-7)(16pi³/15)Gamma7/c^7
           =6/(pi^4 c^7).                                  (7)

Therefore a second complete I bound is

    I<=B_sharp=6^(1/6)/[840pi^(7/6)c^(7/6)],
       c=2^(11/2).                                         (8)

The primary scalar pilot reports B_sharp approximately4.940777144765e-6, also below the prior R006 budget numeral. This is an analytic calibration using a different structural input, not an independent numerical resampling. It supplies no g. To call it wholly independent of all finite geometry would be incorrect: the currently accepted proof of the sharp floor uses the R003 spectral catalog, although it does not use the R010 fan volume.

## 5. Radial normalization: the stronger Rao-Blackwell budget

This refinement is approved by main for the same finite R016 scalar batch; this theory owner still performs no numerical evaluation. Write Y_i=1/(4S_i), R=sum_iY_i~Gamma(32,1), P_i=Y_i/R, and A(P)=sum_i v_i v_i^T/P_i. Then Sigma=A/(4R). With Z_P=detA^(-1/14),

    Z=2sqrt(R)Z_P,
    M=D E_P Z_P, D=2Gamma(65/2)/Gamma32,
    EZ^7=128 Gamma(71/2)/Gamma32 *E_P Z_P^7.                 (9)

In particular D^7 is not the radial factor in EZ^7. Equation(1) becomes

    p1(0)=Gamma(71/2)/[pi^(7/2)Gamma32] E_P Z_P^7.

Hölder on the probability space P alone is also legal and gives a slightly stronger M upper bound. After inserting it into(3),

    I<=eta B_origin,
    eta=Gamma(65/2)/Gamma32
                /[((65/2)(67/2)(69/2))^(1/6)] <1.           (10)

To check both its algebra and direction, put a=E sqrt(R)=Gamma(65/2)/Gamma32 and b=(65/2)(67/2)(69/2). Then ER^(7/2)=ab. The ratio of the P-only Hölder M bound to the raw bound in(2) is a/(ab)^(1/7); raising to7/6 gives exactly a/b^(1/6)=eta. Hölder is strict for the nonconstant R, so a^7<ER^(7/2)=ab, hence a^6<b and eta<1. All moments are finite. No independence other than the established R/P independence has been assumed.

Define B_RB=eta B_fan as the main refined budget and retain B_fan as the raw-origin calibration. Both are strict theoretical upper-budget expressions, with B_RB<B_fan; neither is a lower bound for actualI. This improvement was not numerically evaluated by this theory owner. The raw-origin formula(6) remains independently valid and must not be overwritten or attributed to a P-only calculation.

## 6. Conditional cost, thresholds and what remains unproved

For any separately certified legal full cube/face lower floor g>0, the original positive-source comparison gives

    C7<=7^(7/6)B/g^(7/6),                                  (11)

where B may be the new proved complete I upper budget. The original generalL1, all positive common radii and centered-cube hypotheses are unchanged. The finite-source, obstacle and positive-part mass arguments remain the accepted R001 ones; this round replaces only I.

For a fixed selected budget B the expression in(11) is <=128 exactly when

    g>=g_req(B)=7(B/128)^(6/7).                             (12)

The primary raw-fan pilot gives g_req approximately2.65306467e-6. If a future proof established g=1/320000, the primary pilot's conditional raw-fan expression would be approximately105.744. Neither number is a proved global g or a new C7 bound. A rigorous numeric certificate should use an upper rational B_U in both(11)-(12); the exact test for the conditional expression at g=1/320000 to be below128 is

    7^7 B_U^6 *320000^7 <128^6.

R015's true point upper bound still excludes its larger floor targets. Its fixed old numerator obstruction does not apply to this smaller new numerator. Conversely a single true point value above g_req does not establish the needed global floor. This leaves a concrete original-model g problem, rather than a completed operator upper bound.
