# R017: Gaussian uniform-floor proxies and their paired-budget barrier

2026-10-27. Pure Sol theory and read-only source review; no scalar, cache functional, sample, kernel, MC or MATLAB calculation by this owner. The original7D body64/unit-weight/common-time model and the accepted all-scale/generalL1 comparison remain fixed. Terra's already executed8192old-cache new-functional diagnostic is separate numerical evidence, not a rigorous expectation certificate.

## 1. Definitions and the two distinct lower-bound arguments

For P~Dirichlet(1/2)^64 write A=sum_v vv^T/P_v, B=A^-1, d=(detB)^(1/2),

    C=1/(2pi^(7/2)), D=2Gamma(65/2)/Gamma32, KD=C D,
    q=7/6, F=1/[105*2^(25/6)*pi^(13/12)].

The accepted exact radial representation of the true full-time potential is

    G(x)=KD E_P[d(x^T Bx)^(-3)],                              (1)

in the extended positive sense. Define

    L_vertex(P)=KD d/[max_(sigma in {+/-1}^7) sigma^T Bsigma]^3,
    L_trace(P)=KD d/(trB)^3.                                 (2)

For fixed SPD B, convexity in each coordinate puts the maximum of x^T Bx over the closed cube[-1,1]^7 at a vertex. Therefore L_vertex is a true per-shape cube-uniform lower contribution of the kernel in(1). Sign averaging of vertex quadratics gives E_sigma sigma^T Bsigma=trB. Hence max_sigma>=trB and

    0<L_vertex(P)<=L_trace(P).                               (3)

The larger L_trace is generally not a lower bound for the unaveraged single-shape kernel at every cube point. Its uniform lower-bound status requires symmetry before Jensen, as follows.

The law of P is unchanged when a signed permutation permutes the64unoriented body lines, so the law of B is conjugation-invariant under B7 and d is unchanged. Thus(1) equals the expectation of the finite group average. For x!=0, Jensen for the convex decreasing function u^-3 gives

    avg_Q (x^T Q BQ^T x)^(-3)
      >=[avg_Q x^T Q BQ^T x]^(-3)
       =[|x|² trB/7]^(-3)>=(trB)^(-3), x in[-1,1]^7.

At x0 the original potential is infinite and the lower inequality is immediate. Positive expectations justify averaging even before finiteness is known. Therefore

    G(x)>=g_trace:=E_P L_trace(P)>=g_vertex:=E_P L_vertex(P)
                  for every cube point.                     (4)

Sign averaging alone also suffices for the final cube inequality, since its mean quadratic is sum_j x_j² B_jj<=trB. The full group calculation above gives the more symmetric radial expression; neither approach asserts L_trace is an unaveraged shape floor.

## 2. A deterministic envelope for the diagnostics

By weighted Cauchy-Schwarz, for any u in R7 and positive P with sumP1,

    u^T Au=sum_v (v·u)²/P_v >=[sum_v|v·u|]²=psi(u)².

The accepted sharp floor psi(u)>=sqrt2048|u| therefore proves A>=2048I, B<=I/2048. Determinant AM-GM gives trB>=7(detB)^(1/7). Consequently

    L_trace<=KD (detB)^(1/14)/343
           <=KD/[343sqrt2048].                              (5)

Together with(3), this is the same deterministic envelope for both proxies at every shape, not just the cached8192shapes. It makes their expectations finite. Sampling a bounded quantity does not itself supply a deterministic or probabilistic expectation certificate without its own error argument.

## 3. The exact moment upper bound on either proxy expectation

Let M=E det(Sigma)^(-1/14) be the original raw Gaussian-mixture determinant moment. The established R/P radial identity is

    M=D E_P (detB)^(1/14).

Taking expectation in the first inequality of(5) yields

    g_vertex<=g_trace<=C M/343.                              (6)

This is an upper bound for these specific lower-surrogate values. It is not an upper bound on the actual optimal cube floor infG. Equation(4) and equation(6) must not be reversed or combined into an upper bound on G.

Any pointwise cube-constant minorant chosen before averaging the raw single-shape kernel is at most that shape's exact minimum L_vertex, so its expectation also obeys(6). The present argument also covers the named trace-Jensen proxy. It does not exclude a stronger construction that keeps x-dependence through expectation or takes the exact infimum of a richer symmetry-averaged shape kernel; such a quantity need not be bounded by L_trace.

## 4. A universal barrier for this paired certificate chain

The accepted whole-boundary moment argument gives I<=F M^q. Suppose it is implemented with any certified M0>=M and chooses the corresponding numerator I_cert=F M0^q. Suppose the chosen lower-floor certificate satisfies 0<g_cert<=g_trace, or belongs to a smaller proxy in section3. Then the original positive-source cost expression is

    Cost_cert=7^q F M0^q/g_cert^q.

Using(6) and M0>=M>0,

    Cost_cert>=7^q F(343/C)^q (M0/M)^q
              >=K_barrier:=7^q F(343/C)^q
                =7^(14/3)pi³/840.                           (7)

The constants can also be checked without fractional powers of pi by using F=C^q H/896, H=16pi³/15: then K=(7*343)^q H/896=7^(14/3)pi³/840. Thus there is no dependence on how tightly M is estimated.

The strict inequality K>128 has a simple exact proof requiring no scalar calculation: 7^(2/3)>3 because49>27, and pi>3, so

    K>7^4*3*27/840=194481/840>128,

the last comparison being194481>107520. The elementary pi>3 follows, for example, from the inscribed regular hexagon in the unit circle. The later four-scalar batch is a calibration and useful scale report, not a missing numerical premise of this exclusion.

This is a barrier for a selected proof chain and its reported certificate expression. I_cert is an upper budget, not an identity for actualI. A lower bound on this expression is not a lower bound for C7. It does not refute the original Green model, an x-dependent method, a different numerator inequality, or all conceivable uniform Gaussian bounds. Even at an exact M, the named proxies combined with the same Minkowski numerator cannot certify128.

## 5. The R016 finite scalar specialization

R016 proved a sharper M upper bound by origin moments and Rao-Blackwellization and, through the same F M0^q interface, produced the rational I budget

    J0=4018421141119/10^18.

Define M0=(J0/F)^(1/q). This is a legal moment upper bound here because the R016 underlying expression was F times a proved moment upper bound to powerq and J0 rounds that expression upward. Merely knowing an arbitrary unrelatedIupper would not establish this fact.

The named proxy expectations are therefore at most U_proxy=C M0/343, whereas the same J0 budget needs g_req=7(J0/128)^(1/q) to reach128. Their ratio is exactly

    U_proxy/g_req=(128/K_barrier)^(1/q)<1.                   (8)

Hence an improved finite-cache estimate of the named proxies cannot repair this same-budget route. Terra's already executed means L_vertex about1.935e-8 and L_trace about1.439e-7 are descriptive new functionals of the old cache; they are not used to prove(7)-(8), and no new probes are needed.

## 6. Provenance and scope

The next authorized batch is only four deterministic512-bit scalar enclosures specified in registration.md. It may use the saved R016 rational budget, C,D,F and the sharp floor. No cache reread for new functionals, MC, point search, Green grid, optimizer, old baseline, MATLAB or full kernel is authorized. The current theory owner computes no scalar. Parent retains proof responsibility and integration.
