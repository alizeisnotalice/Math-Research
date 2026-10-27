# R015: true Green point, local continuity and the scope of a failed fixed budget

2026-10-27. Sol theory and read-only implementation review. No numerical kernel/phase, threshold arithmetic, MATLAB or MC run by this owner. The actual numerical enclosure/comparisons belong to the separately executed green_point certificate; its sources are hash locked. The model is still the original7D body64, unit weights and common Cauchy time. The accepted complete R010 fan and R011 Fourier/DD identities are assumed with their existing structural certificates.

## 1. Exact point and the finite nonzero hypothesis

The physical point x0 is the exact dyadic conversion of R003 Sobol index2710, not a decimal approximation:

    0x1.0000000000000p+0
    0x1.cb8b8cf800000p-1
    0x1.ea76a78800000p-1
    0x1.d94dfdf800000p-1
    0x1.18baa38000000p-5
    0x1.7babca4800000p-1
    0x1.d943b29000000p-1.

Thus x0_1=1, all other coordinates are strictly between0 and1, and all seven absolute coordinates are distinct. These assertions follow by exact binary values and are also checked by the executed numerical certificate. If x lies in a span of k<=3 body directions, each coordinate is a signed sum of the same k coefficients. There are at most2^(k-1)<=4possible absolute values, since antipodal sign patterns have the same absolute value. Consequently this point is off S3. This is not a new classification of S3 or a rerun of that catalog.

Let u_r=xi_r/q_r be the46fixed normalized fan rays, and b_(Q,r)(x)=x·Q u_r. The source R012 whole half-orbit certificate asserts every one of its322560*46exact ray projections is nonzero. Opposite Q changes the sign, so this is the same nonzero hypothesis for the complete B7group. R015 independently guards all46integer projections at every executed phase as well. This is a finite exact test, not an inference from floating magnitudes or a measure-one genericity claim.

## 2. Epsilon0 is the true G, not a small-cutoff surrogate

For a top simplex S with alpha_S=|detW|/product q and for epsilon>0, the accepted complete Fourier formula is

    H_epsilon(x)=(2pi)^(-7) sum_Q sum_S alpha_S
                 integral_epsilon^infinity product_(r in S)
                                          (t-i b_(Q,r)(x))^(-1)dt.

It equals integral_epsilon^infinity p_t(x)dt for the actual nonnegative Cauchy density. At x0 all b are nonzero. Every individual cone time integrand is absolutely integrable on(0,infinity): bounded near0 and O(t^-7) at infinity. The group and simplex sets are finite. Dominated convergence on these cone integrals and monotone convergence of the actual positive density therefore give the same limit,

    G(x0)=integral_0^infinity p_t(x0)dt
         =2/(2pi)^7 sum_(Q in half B7) sum_S alpha_S
                        DD6[log|b|](b_(Q,r)(x0):r in S).        (1)

The group stabilizer at x0 is trivial because all absolute coordinates are nonzero and distinct. The factor2 is opposite-pair conjugacy; it is not an assumed symmetry of an individual ordered simplex. Equation(1) has no shape-expectation, frequency-truncation or omitted short-time error.

For sorted exact real nonzero knots, DD uses the ordinary difference recurrence and, for a repeated subblock at b, the Hermite coefficient

    (log|b|)^(k)/k!=(-1)^(k-1)/(k b^k), k>=1.                 (2)

Mixed positive and negative knots are permitted. Their distinct-knot differences stay exact; coincident-knot derivatives are taken locally away from0. The formula follows equally by taking real parts of the right-half-plane epsilon>0 logarithms and then epsilon down to0. No analytic branch through b0 is assumed.

With x_j=n_j/d_x, L=lcm q_r, B=d_xL, the exact knot is a_r/B for integer a_r. Since log|a/B|=log|a|-logB and a positive-order DD annihilates constants,

    DD_b^6 log|b|=B^6 DD_a^6 log|a|.                           (3)

The program must multiply by B^6 once, not at every recurrence order. A zero a is a failure of this sufficient formula; it must never be deleted, replaced by an epsilon, or classified by itself as a physical Green singularity.

## 3. A full ambient neighborhood of continuity

Continuity can be proved directly here without appealing to R001's stronger global off-S3 theorem. Finiteness of Q,r and their exact nonzero values implies

    delta=min_(Q,r)|b_(Q,r)(x0)|>0,
    R=max_r |u_r|>0.

For |x-x0|<delta/(2R), every |b_(Q,r)(x)|>=delta/2, since signed permutations preserve Euclidean lengths. For0<t<=1,

    |product_(r in S)(t-i b_(Q,r)(x))^(-1)|<=(2/delta)^7;

for t>=1 the same absolute product is at most t^-7. This common bound is integrable over all t>0 and uniform over this ambient neighborhood and every cone. Each integrand is continuous in x. Therefore every cone integral, their finite sum and hence the actual G are continuous throughout that open neighborhood. The same hypotheses show G is finite there.

This argument establishes a positive neighborhood radius abstractly; the point calculation did not compute a quantitative spatial cover. No numerical minimum-delta estimate or new kernel point is needed to deduce continuity or the essential-floor consequence below.

## 4. Boundary points constrain essential floors

Let Qcube=[-1,1]^7. Suppose g is a legal essential cube floor: G(x)>=g for Lebesgue-almost every x in Qcube. Then

    g<=G(x0).                                                 (4)

Otherwise continuity provides a ball around x0 on which G<g. Its intersection with the interior of Qcube has positive7Dvolume: move the first coordinate slightly below1 and take a sufficiently small interior ball; the other six coordinates already lie strictly between0 and1. This contradicts the a.e. lower bound. Thus the fact x0_1=1 does not create a measure-zero escape.

For a legal essential floor on the positive face F={(1,y):y in[0,1]^6}, the same conclusion holds. The other six coordinates of x0 are strictly interior, so the continuous restriction has a relative-open neighborhood of positive6Dface measure on which G<g if g>G(x0). A literal pointwise full-face/cube floor is of course also constrained. Signed-permutation invariance then transfers the analogous statements to other congruent faces, but is not needed for this contradiction.

## 5. What the executed outward endpoint excludes

The separate numerical certificate records an outward dyadic interval [L_G,U_G] for the true G(x0). Its displayed rigorous decimal brackets are

    0.000003782255829017294943880908 < G(x0)
                       <0.000003782255829017294943880909.

The exact upper endpoint used for acceptance is U_G=N/D in green_point/assessment.json, not a decimal midpoint. The numerical owner actually ran the exact Fraction tests U_G<1/K for K219000,194000,179000 and the accounting checker. This theory owner read their saved results and source code, but did not rerun those arithmetic tests or any phase. Acceptance reduces to the transparent integer comparisons K*N<D. As long as the locked outward certificate/checker is accepted, (4) proves all three proposed essential/global floor targets impossible for this G. It does not give any operator C7 lower bound or refute a smaller floor target.

For reproducibility the locked endpoint is

    N=3587729854234050249868222360880233704712469649049936819188541798427025246713,
    D=948568795032094272909893509191171341133987714380927500611236528192824358010355712.

No old H interval is used as a G value. The strict consistency H_(1/1024),upper<L_G, checked by the numerical owner, is useful but does not prove(1); the cone limit/MCT argument does.

## 6. The fixed numerator budget and its limited exclusion

Fix the specific already proved complete-I upper budget, as a chosen numeral,

    I0=7723627498344/10^18.

For a legal positive g the associated certificate expression is

    C_expr(g)=7^(7/6)I0/g^(7/6).

By(4), g<=G(x0)<=U_G, and hence

    C_expr(g)>=7^(7/6)I0/U_G^(7/6).                            (5)

Its sixth power is the exact rational 7^7 I0^6/U_G^7. With I0=A/B, the test that the right side of(5) exceeds128 is

    7^7 A^6 D^7 >128^6 B^6 N^7.                              (6)

The separately executed assessment records this strict integer/Fraction comparison as true; the displayed sixth-root value is approximately158.38947. The root display is not needed for the exact acceptance test.

Equation(5) is a lower bound on this fixed certificate expression, not on the actual optimal operator constant. I0 is an upper budget for I, so I<=I0; a smaller true I or a sharper proved upper budget can reduce the numerator. In particular this argument does not prove actualI is too large, does not exclude the entire original Green model, and does not prove C7>=128 or C7>=158.38947. Its concrete conclusion is that combining this fixed numerator I0 with any legal g from this same original G cannot certify a bound at most128. Progress using this route must improve the numerator or change the mechanism, not pursue one of the already excluded floor targets.

## 7. Read-only implementation audit

The actual path is green_point/run.py importing green_kernel.py. Compared with the accepted R012 kernel, the live generic path changes only the logabs/repeated-derivative branch, asserts epsilon0 plus all46exact nonzero projections, and calls the generic phase with epsilon0. The real DD recurrence, common denominator, B^6 factor and2/(2pi)^7 phase normalization are unchanged. The repeated-knot probe at integer1 correctly has DD6logabs=-1/6.

All5040Gchunks are new. Pilot chunks0..3 are computed as G, pass a planning-width gate, and remain in the same final sum exactly once. Remaining chunks4..5039 have checked phase IDs and counts. The final sum uses the already normalized signed chunk intervals and applies no second factor. Exact source/hash checks and per-phase zero-node guards are present; no clipping or replacement appears. Assessment/checker compare the saved outward dyadic endpoint, not a midpoint, and use the sixth-power inequality in the correct direction.

green_kernel.py retains unused old H scalar/origin/main helper functions. They are not valid epsilon0 entry points, and the actual run.py does not call them. Their presence is not evidence of new G scalar/origin tests. This is a documentation/reproduction limitation, not a blocker for the executed generic path. No prior sealed files were changed by this reviewer.
