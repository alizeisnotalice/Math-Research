# A008 alpha=3/2: a strict global cube g lower bound, still noncompetitive

For the existing equal-weight alpha=3/2 body32 model, the batch proves

    ess inf_(x in [-1,1]^6) G_(3/2)(x)
       > 54117529919/500000000000000000
       = 1.08235059838e-7.

The proof covers the entire cube, including its boundaries and singular
strata as an extended-valued pointwise lower inequality. It uses no alpha=1
kernel values or Monte Carlo minimum. The best of the four registered
analytic schemes is a global degree-6 cosine lower polynomial, with an
improved positive zeroth-order integral from spherical Jensen.

The bound is far too small to combine competitively with the current F008
whole-face norm budget. All four optimized allocated costs exceed 1200000;
the best is about 1235395.39. This is a failure of these allocated budgets,
not a lower bound on actual operator/model cost or an upper bound on actual
g. C6<=64 remains the applicable upper endpoint.

## Exact second-order construction

With the Fourier normalization (2pi)^-6 and

    c6=32(8/3)^(-1/4),
    c6 ||xi||^(3/2) <= psi(xi) <=32 ||xi||^(3/2),

the global inequality cos(s)>=1-s^2/2 gives

    p_t(x) >= a t^-4-b ||x||^2 t^(-16/3),
    a=Gamma(4)/(96 pi^3 32^4),
    b=Gamma(16/3)/(1152 pi^3 c6^(16/3)).

The two signed terms use opposite spectral envelopes. There is no claim
that Fourier ordering alone implies a pointwise density comparison.
On the cube, the correct Euclidean bound is ||x||^2<=6. Thus

    G(x) >= a/(3T^3)-(18b/13)T^(-13/3).

Its unique optimal positive cutoff is T=(6b/a)^(3/4), and the resulting
global lower certificate is

    g_2=(4a/39)*(a/(6b))^(9/4).

Arb384 gives an enclosure around 2.252927563614421145e-8; the cutoff is around
0.20607542962918305. Both the closed formula and an independent polynomial
root optimization agree.

## Improved zeroth-order term and two finite higher orders

The exact spherical mean of psi on S5 is

    mu=64*6^(3/4) Gamma(5/4)/(sqrt(pi) Gamma(15/4))
       =28.367675982807811148...,

whereas c6=25.041353281172597065... and the pointwise upper rate is 32.
Gaussian/polar integration proves the mean formula. Angular Jensen yields

    average_omega exp[-t r^(3/2)psi(omega)]
      >=exp[-mu t r^(3/2)].

Consequently the positive zeroth-order coefficient may use mu instead of
32. This is an angular integral estimate; psi<=mu is not asserted.

For global cosine lower polynomials of degree 2,6,10, let K=1,3,5 and

    A_k=Gamma(4+4k/3)
        /[96 pi^3 4^k k! (3)_k d_k^(4+4k/3)].

Use d_0=mu in the improved schemes, d_k=32 for positive even k>=2, and
d_k=c6 for negative odd k. Then for every x and every t>0,

    p_t(x)>=sum_(k=0)^K (-1)^k A_k ||x||^(2k)t^(-4-4k/3).

The cosine polynomial bounds are global, proved by successive second
derivatives and integration from zero, so there is no omitted frequency
region or unbounded Taylor remainder. All Fourier and positive-time tail
integrals converge absolutely where they are used.

For r=||x||>0 take T=r^(3/2)z^(-3/4). Positive integration gives

    G(x)>=r^(-9/2)F_K(z),
    F_K(z)=z^(9/4)sum_(k=0)^K (-1)^k A_k z^k/(3+4k/3).

For positive F_K, replacing r by sqrt(6) gives the cube lower bound
6^(-9/4)F_K(z). At x=0, G(0)=+infinity, so it causes no exception to the
lower inequality. The self-contained proof in `theorem.md` also obtains
p_t>=0 by truncated compound-Poisson measures, rather than importing a
Cauchy density comparison.

The derivative is exactly

    F_K'(z)=(3/4)z^(5/4)sum_(k=0)^K (-1)^k A_k z^k.

For degrees 6 and 10, the recorded strict coefficient inequalities make the
last polynomial strictly decreasing for z>=0. Each therefore has one positive
root, giving the global optimum of this one-variable certificate. Arb384
strictly brackets each root with 220 rational bisections. A rational point
in the bracket independently supplies the constructive g certificate.

## Actual finite results

|Scheme|Optimized whole-cube certificate value|Optimal cube tail cutoff|Allocated cost with fixed F norm budget|
|---|---|---|---|
|Degree 2, radial 32/c6|2.2529275636144211e-8|0.20607542963|10014600.337|
|Degree 2, spherical Jensen|1.0789171484418263e-7|0.14356466181|1240640.056|
|Degree 6, spherical Jensen|1.0823505983878280e-7|0.14487556321|1235395.394|
|Degree 10, spherical Jensen|1.0818122615104217e-7|0.14494865046|1236215.147|

The complete Gamma coefficients, strict root brackets, constructive rational
z values, and Arb enclosures are in `results.json`. `checked_results.json`
also provides narrow strict rational intervals for the optimized certificate
values. In particular the best lies in

    (54117529919/500000000000000000,
     108235059839/1000000000000000000).

This interval describes the optimized lower-certificate formula, not the
unknown true essential infimum. Increasing the cosine degree did not improve
the fourth scheme; no further order was tried.

For the degree-6 method, the exact coefficient formula has Arb values near

    A0=3.11268640077026263718e-9,
    A1=3.89984656777506942396e-11,
    A2=3.14487456662392713434e-14,
    A3=4.75296754252996918936e-16.

The strict monotonicity check A2^2/(3 A1 A3)<1 has value about 0.0177858.
The degree-10 checks have ratios about 0.0533574 and 0.00490250, also strictly
below one. No stationary root was dropped by a numerical optimizer.

## Budget capability and its exact scope

F008's actual norm-envelope formula has value

    U=0.001002680234066613134146...,

and the existing general-L1/all-positive-scale alpha=3/2 interface has
allocated cost (4U/g)^(4/3). To make that expression less than 64 requires

    g>U/(4sqrt(2))
      =0.000177250498217554217072277498... .

The independent reconstruction verifies every one of these four optimized
certificate values is below 109/10^9, while the required threshold is above
177/10^6. Equivalently every optimized allocated cost is greater than
1200000. Hence no choice of tail cutoff within any of these four schemes
can close the cost with this fixed U budget.

This conclusion concerns the formulas used to allocate upper and lower
certificates. The actual face norm can be much smaller than U, and the actual
g can be much larger than the proved lower bound. Neither is bounded in the
opposite direction here. No old MC minimum is promoted to an exact statement.

## Execution, verification, and retained scope

Exactly four analytic schemes and their four cost-budget scenarios were run,
matching the eight-scenario registration. `check_bounds.py` completed with
exit0 in **0.089635333 seconds** of script time (tool chunk b9cddb); its actual
Python PID is in the results. There were zero kernel evaluations, zero MC
samples, no 5D grid, no MATLAB calls, no subagents and no pushes.

SymPy actually verified the tail derivative identities and the cosine
second-derivative induction identities through degree10. A separate
`verify_normalization.py` reconstructs all angular coefficients using the
alternative Gamma moment formula, checks the base Fourier normalization,
reconstructs root signs and rational lower bounds, and proves the stated
budget comparisons. It introduces no new parameter scenario.

C008 coercivity and F008 norm/interface files were read-only and SHA-256
locked. After computation, the proof's positivity explanation was expanded
without changing any coefficient or value; both the executed proof hash and
final proof hash are recorded. The exact alpha=3/2 symbols are used throughout.
The completed chamber counterexample and round007 files were untouched.

Reproduce with:

    /Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python check_bounds.py
    /Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python verify_normalization.py

The batch is complete and stopped. Retained artifacts are the finite
registration, self-contained theorem, this report, the two scripts, input
manifest, full results, and independent reconstruction. No C6 endpoint change
is registered by this subtask.
