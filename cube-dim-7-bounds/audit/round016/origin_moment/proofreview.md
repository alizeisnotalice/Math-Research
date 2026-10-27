# R016 independent origin-moment and radial-conditioning audit

2026-10-27. **Accepted.** The original complete integral has the strict rational upper budget

    I <4018421141119/10^18.

This is the radial-conditioned origin-moment budget, not a new global floor or operator bound. The independent proof check and the registered768bit ten-scalar audit found no mathematical/source-interface blocker. The parent retains overall proof and integration responsibility. Earlier sealed rounds were not changed.

## 1. Scope, provenance and assumptions

This audit used a new no-hardlinks partial/sparse clone at14f3507f406185667f84146328f1db9da3f9949a and wrote only audit/round016/origin_moment. Registration preceded arithmetic. The accepted original model has64body columns v_i in R7, unit weights and one common Cauchy time; the columns span R7. The R001 original-G Gaussian/boundary comparison and generalL1/allr source argument, accepted R009/R010 full fan, and accepted R011 sharp-floor certificate are explicit dependencies. No new fan, spectral catalog, phase, DD, shape sample or spatial cover was executed. Reading a prior source/hash is not claimed as a new independent proof of its enumerative completeness.

The current Sol theory proof was fully read, alongside the actual numerical run.py/check.py, manifest, pre-run hashes, receipt, result and report. The independent checker imports no author module and runs no author checker. It validates the new numerical manifest and all members; theory manifest and all members; all registered pre-run input hashes; copied R010cells and R006certificate against their copied accepted manifests; and sharp-floor/primary-pilot records against their source hashes. The spectral catalog and R011pilot result are hashed, not mathematically re-enumerated or re-evaluated. Exactresults.json contains24source hashes; sources.json preserves them separately.

The prior primary scalar pilot, Terra's actual512bit ten-scalar run and this independent768bit ten-scalar run have separate ownership. This audit's sole checker invocation completed exit0 on its first attempt, with empty stderr. No failure, repair or scientific rerun occurred.

## 2. Pointwise origin identity and the raw seventh moment

Let S_i>0 be independent with Laplace transform E exp(-lambda S_i)=exp(-sqrt(lambda)), and Sigma=sum_i S_i v_i v_i^T. Since the columns span, Sigma is SPD a.s. For t>0,

    E exp(-t^2 xi^T Sigma xi)=exp(-t sum_i|v_i dot xi|)=exp(-t psi(xi)).

This preserves the common-time generator. The conditional Gaussian covariance is2t^2Sigma, whose density is

    (4pi t^2)^(-7/2) detSigma^(-1/2)
            exp[-x^T Sigma^-1 x/(4t^2)].

At the actual selected point0 the Fourier integrands are nonnegative. Positive Tonelli and the Gaussian integral give

    p1(0)=(2pi)^(-7) integral exp(-psi(xi)) dxi
         =E[(2pi)^(-7) integral exp(-xi^T Sigma xi) dxi]
         =(4pi)^(-7/2) E detSigma^(-1/2).

The original characteristic function is integrable (a spanning finite family gives psi>=c|xi| with some c>0), so its inverse Fourier density is a continuous canonical version. This identity establishes its pointwise value directly; equality of density versions only a.e. is not used at the selected point. The expectation is finite by the same frequency bound. No finite expectation of Sigma is needed. The sharp floor later gives a quantitative bound, but is not necessary just for finiteness.

Put rawZ=detSigma^(-1/14) and M=E rawZ. Then rawZ^7=detSigma^(-1/2), so Hölder on the probability space gives

    M<=(E rawZ^7)^(1/7)=2sqrt(pi) p1(0)^(1/7).

The exponent and direction are correct. p1(0) is a finite time-one density, whereas p_t(0)=t^-7p1(0) and G(0)=infinity. The latter does not invalidate the time-one identity or the finite complete-face integral obtained below.

## 3. Full-boundary normalization, including every power of2 and pi

Independently integrating the conditional density over all positive time yields, for Q=x^T Sigma^-1x,

    integral_0^infinity t^-7 exp[-Q/(4t^2)]dt=64Q^-3,
    G_Sigma(x)=detSigma^-1/2 Q^-3/(2pi^(7/2)).

The factor64combines with(4pi)^(-7/2)=1/(128pi^(7/2)), leaving the required1/(2pi^(7/2)). For q=7/6 the entire boundary of[-1,1]^7 satisfies

    integral_boundary (x^T Sigma^-1x)^(-7/2)dσ
        =|S6| sqrt(detSigma),  |S6|=16pi^3/15.

This follows by radial coordinates on the cube: x dot n=1 on every face, and integrating the homogeneous degree-7 angular function maps surface area to the spherical one. Applying L^q Minkowski to the full boundary first gives

    ||G||_(Lq(boundary)) <= (2pi^(7/2))^-1 (16pi^3/15)^(6/7) M.

The determinant exponent is-1/2+3/7=-1/14, exactly rawZ. Only the averaged original G is invariant under coordinate sign/permutation changes. Its boundary q-integral is896times I:14faces times64orthants per face. There is no coordinate symmetry assumption on a fixed conditional Gaussian. Raising the preceding norm inequality to q and dividing by896 gives

    I <= M^(7/6)/(105*2^(25/6)*pi^(13/12)).

Positive truncations followed by monotone convergence justify the inequality before finiteness has been concluded. Substituting the raw moment bound cancels powers exactly and yields

    I <= B_raw=p1(0)^(1/6)/(840sqrt(pi)).

This is a complete-time, whole-face original-G budget. Only the numerator estimate changes in the accepted positive-source comparison; the generalL1, all positive common radii, centered-cube and finite-source/obstacle passage retain the R001 proof. No new truncated-kernel comparison or generator is being introduced.

## 4. Exact fan density and the independent sharp-floor calibration

In each top cone, psi(sum_i a_i xi_i)=sum_i a_i q_i for a_i>=0. The determinant Jacobian and seven elementary positive exponential integrals give

    integral_cone exp(-psi(xi)) dxi=|detW|/product_i q_i=alpha.

Here xi=Mw with detM=1. This is an exponential cone integral, not a unit-sublevel simplex volume; no additional7! or q-normalization belongs in alpha. The accepted fan supplies212positive alpha values over the ordered chamber, and the645120signed permutations cover allfrequency space up to measure-zero boundaries. Consequently

    p1(0)=645120 sum(alpha)/(2pi)^7=5040 A_sum/pi^7,
    A_sum=907157967345860903220972676779187
          /28127466839029475035590994122592223232000000000.

The independent Fraction ledger verifies135cells,212distinct7vertex simplex keys, positive q and alpha, and each saved alpha=abs(detW)/product(q), then sums them exactly to the stated rational. It does not recompute geometric determinants or pulling coverage; those are accepted prior structural certificates. The copied file hash binds this exact sum to that certified fan.

Thus B_raw=(5040 A_sum)^(1/6)/(840pi^(5/3)). Independently of this fan-volume functional, the accepted sharp floor c=32sqrt2 and radial integration give

    p1(0)<=|S6|Gamma7/[(2pi)^7 c^7]=6/(pi^4 c^7),
    I<=B_sharp=[6/(pi^4 c^7)]^(1/6)/(840sqrt(pi)).

Its current proof uses the accepted spectral catalog and previously executed21orbit maximum-norm receipt. It is independent of the R010volume sum, not independent of all finite geometry. The new scalar intervals separate p_fan strictly below this calibration. Both budgets are upper budgets; a lower endpoint of either budget expression is not a lower bound for actualI.

## 5. Radial conditioning: strict eta and the correct seventh factor

Use Y_i=1/(4S_i), iid Gamma(1/2,rate1). Their sum R is Gamma(32,rate1), independent of P_i=Y_i/R. Let A(P)=sum_i v_i v_i^T/P_i and Z_P=detA(P)^(-1/14). Then

    Sigma=A(P)/(4R), rawZ=2sqrt(R) Z_P,
    M=2a E_P Z_P, a=Gamma(65/2)/Gamma32,
    E rawZ^7=128 E R^(7/2) E_P Z_P^7
             =128 Gamma(71/2)/Gamma32 E_P Z_P^7.

The radial factor in the last line is not[2a]^7. This distinction is material; radial Hölder loss can be removed before applying Hölder to P. Let b=(65/2)(67/2)(69/2). Gamma recurrence gives E R^(7/2)=ab. Hölder applied only to P improves the raw M estimate by factor a/(ab)^(1/7), and raising to7/6 in the I interface gives

    I <= B_RB=eta B_raw,
    eta=a/b^(1/6).

Strict Hölder for the nonconstant positive R gives a^7<E R^(7/2)=ab, hence eta<1. All moments in use are finite. This strict improvement does not depend on Z_P being nonconstant. Finiteness of E rawZ^7 and independence imply the needed finite P moment. The process itself and its common time remain unchanged.

## 6. Independent768bit scalar path and exact acceptance

The auditor uses an exact half-integer reduction, different from the producer's direct Gamma and log/exp path:

    a=C sqrt(pi), C=64!/(4^32 32!31!),
    Gamma(71/2)/Gamma32=C b sqrt(pi).

The latter recurrence is verified using exact factorial fractions. For the main bounds the independent positive-root expressions are

    B_raw^6=5040 A_sum/(840^6 pi^10),
    eta^6=C^6 pi^3/b,
    B_RB^6=5040 A_sum C^6/(840^6 b pi^7).

Only ten named scalars are computed at768bits: p_fan, sharp p upper, raw M upper, B_raw, B_sharp, eta, B_RB, the required-g expression from rounded B_RBupper, and the two conditional costs from rounded rational budgets. There are8explicit positive-root operator calls,2square-root calls and1pi call in the scalar path; log/exp/Gamma API calls are0. These are ten outputs, not ten elementary/transcendental operations.

All ten independent768bit intervals lie within their corresponding frozen512bit intervals. This is stronger than overlap consistency, while validity still rests on the proved expressions and directed arithmetic. The primary pilot's raw I interval also contains the independent raw expression. The exact outward endpoints are saved as full rationals in exactresults.json. Representative values are

    p1(0)≈5.3818764733219049339e-14,
    B_raw≈4.1270024988294296110e-6,
    B_sharp≈4.9407771447653848703e-6,
    eta≈0.97369001890802691650,
    B_RB≈4.0184211411187016502e-6.

Ceiling the independent upper endpoints at denominator10^18 reproduces the exact producer budgets

    B_raw <4127002498830/10^18,
    B_RB <4018421141119/10^18.

These rounding choices are fixed arithmetic, not optimization. The upper endpoints are strictly below the grids. The RB interval lies strictly below the raw interval; the raw rational upper lies strictly below the saved lower endpoint of the old R006 deterministic upper-budget expression. B_sharp<4940777144766/10^18 is also below that old expression. The direction compares two proved upper-budget expressions; it never infers a lower bound on actualI, actualM or C7.

## 7. Conditional floor and remaining gap

For a separately proved legal g>0 the accepted comparison would give

    C7 <=7^(7/6) B/g^(7/6).

Using the rounded B_RBupper, the required floor for the expression to be<=128 is g_req=7(B_U/128)^(6/7), with independent value approximately2.59312073253869755e-6. At the hypothetical exact g=1/320000, the independent exact Fraction tests give

    7^7 B_RB,U^6/g^7 <7^7 B_raw,U^6/g^7 <128^6.

Their sixth roots are approximately102.9622769241 and105.7444103613. These are conditional certificate expressions, not new C7bounds: this round supplies no globalg. R015 still excludes the larger g targets it tested. Its fixed OLD numerator cost obstruction does not obstruct the new smaller numerator, and its single true point above this new hypothetical floor proves neither legality nor illegality of that floor.

## 8. Reproduction and stop

minchecker.py is stdout-only and reads explicit CLI paths: --numeric (frozen R016origin_moment), --theory (frozen R016origin_moment_theory), --primary (primary R016scalar JSON), --sharp (prior sharp-floor JSON), --catalog (its spectral-catalog source JSON), --pilot (the primary record's R011results source). It does not use cwd, author generator imports or the absolute paths embedded in those records to find inputs. Supply the paths and redirect stdout externally to save exactresults.json. Existing source basenames are retained for resolving the pre-run list. Python with python-flint0.9.0 is required.

This run used the supplied local venv, exact Fraction and Arb768. Source/metadata checks and scalar arithmetic are distinguished from mathematical proof and prior scientific certificates. There were0phase/DD/MC/MATLAB calls,0new kernel points, no commits/pushes and no changes to sealed sources. No repeated baseline, mesh, optimizer or future task was started. The audit is complete and stops here.
