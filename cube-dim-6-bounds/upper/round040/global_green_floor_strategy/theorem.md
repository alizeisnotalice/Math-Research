# A040: S3 obstruction and a correctly paid finite-domain formulation

This concerns only the original body32 Cauchy model. Write V={(1,s2,…,s6):si=±1}, f_t(z)=t/[π(t²+z²)], p_t=density(ΣvZ_v), G=∫_0^∞p_t dt, and S_j for the union of spans of at most j body directions. The target g*=1/65000 is unproved. No chamber convexity or full-G higher regularity is assumed.

## 1. An explicit S3 stratum has infinite all-time curvature moment

Take the six body columns, in this order,

    A=[1, flip5(1), flip6(1), flip2(1), flip3(1), flip4(1)].

Here flipj changes only coordinate j to −1. Their determinant has absolute value 32. For x=(1,1,1,1,b,c), 0≤c≤b≤1,

    A^−1 x=((b+c)/2,(1−b)/2,(1−c)/2,0,0,0).

The first three columns span U={x1=x2=x3=x4}; it is a genuine S3 plane. Let R contain the remaining 26 body columns. Direct exact structure calculation gives the row absolute sums of A^−1R:

    (49,15,15,15,15,15).

In the positive full-rank coarea identity,

    p_t(x)=1/32 ∫_(R^26) [∏_(i=1)^6 f_t((A^−1(x−Rz))_i)]
                                      [∏_(j=1)^26 f_t(z_j)] dz,

restrict all 26 remaining variables to |z_j|≤t. This event has product mass 2^−26. Write z=t u. For 0<t≤1/100, the first three solved parameters have absolute values at most 149/100,13/20,13/20. The last three are at most 15t. Consequently their Cauchy factors obey

    first three product ≥ t³/[π³ A0 A1²],
    last three product ≥ 1/[π³ t³ 226³],
    A0=(1/100)²+(149/100)²=11101/5000,
    A1=(1/100)²+(13/20)²=2113/5000.

Thus, uniformly over the entire displayed S3 face triangle,

    p_t(x) ≥ c* := 5000³/[2³¹π⁶ 226³ 11101·2113²]
                 >1/10²⁰,                   0<t≤1/100.          (1)

The strict rational comparison was actually checked with Arb384. The formula is a positive change of variables followed by restriction of an ordinary integral, not conditioning on a zero-probability fiber.

It follows for every T>0 that

    ∫_0^T t^−2 p_t(x)dt=∞,  and  ∫_0^T t^−1 p_t(x)dt=∞.         (2)

The moment obstruction is not confined to this example. Given any x∈S3, select three independent body columns whose span contains x, complete them to a six-column body basis, and repeat the argument: the last three solved coefficients at x vanish, the remaining-variable box has fixed positive probability, three factors are ≥const·t and the other three ≥const/t. The constant may depend on the chosen basis and x but is positive. Such a completion exists because the frame spans R6. Any three distinct body lines are independent: relative to a pair's nonempty agreement and difference blocks, a third vector av+bw with all coordinate absolute values one would require |a+b|=|a−b|=1, hence ab=0 and would repeat one of the two lines. Lower-rank cases can first be extended to a triple.

Therefore J2=∞ on all of S3. This statement does **not** prove that an actual directional derivative or Hessian of G diverges: signed derivative cancellations may occur. It proves that the positive coarea bound D_h²p_t≥−|h|²p_t/(16t²) cannot give a finite all-time Hessian constant there. Replacing the positive-band U2 by an all-time J2 at such points is invalid.

For this stratum the failure is quantitative. If J2_(a,T)(x) denotes the positive band moment and 0<a<min(T,1/100), then

    J2_(a,T)(x) ≥ c*[1/a−1/min(T,1/100)].                       (3)

Any valid whole-cell upper enclosure that includes this point must grow at least at this rate as a↓0. Merely shrinking the lower time cutoff cannot preserve a fixed moment-based Taylor radius. This is an obstruction to that certificate route, not to the truth of the desired floor.

## 2. The obstruction survives the actual S2 payment

The frozen A012 formula on the ordered positive face Δ={1=x1≥x2≥…≥x6≥0} is

    d2(x)²=min_(k=1,…,5) q_k(x),
    q_k=Σ_(i≤k)(xi−mean_topk)²+Σ_(i>k)(xi−mean_bottom)².

This is used as a proved input, not re-enumerated. At x=(1,1,1,1,1/2,0), direct exact substitution gives min q_k=1/8. This point is in S3 but is far outside the tiny S2 barrier tube. In particular the S3 difficulty cannot be removed by relabeling it S2.

The proven A009 short-potential lower barrier, combined with G≥K_(4/125), gives throughout the cube

    G(x) ≥ [1/(5·10⁹)][1/d2(x)−40000]_+,

with G=∞ on S2. Set δ0=1/125000. Then the entire **closed** tube d2≤δ0 has the correctly paid floor

    G(x) ≥17/10⁶ >1/65000.                                    (4)

The saved A038 third-point upper endpoint is strictly below 17/10⁶. Hence a cube minimum is not lost by removing this tube; even its boundary cannot be a global minimizer. No factor496 or disjointness assumption occurs in this lower payment.

## 3. A finite semialgebraic domain retains every possible minimum

Define the new target-domain description

    D={x∈Δ : q_k(x)≥δ0² for k=1,…,5}.                         (5)

It is compact, nonempty, and is described by finitely many rational linear/quadratic inequalities. It retains the complete ordered-face boundary, ties and zeros, all S3 strata not paid by (4), and their intersections. It is still five-dimensional. Nothing here reduces the search to first4=1, to the explicit two-dimensional S3 triangle, or to finitely many points.

For completeness, G is finite and continuous off S2 without assuming any higher derivatives. If dist(x,S2)≥δ>0 and Aa=x, at least three parameters exceed Rδ=δ/(60√6) in absolute value; otherwise all but at most two columns contribute at most 30√6 Rδ=δ/2 to the distance from a pair span. Sum the positive coarea submeasures over triples. Each triple is independent as above; its Gram determinant is a positive integer and hence at least one. The remaining-variable quotient is three-dimensional, and its symbol equals ψ restricted to that triple's perpendicular space. The A009 moment proof gives ψ≥8√6|ξ| there. Its density peak is at most t^−3/[π²(8√6)³]. Three selected factors on |a_i|>Rδ are at most t³/(π³Rδ⁶). Thus the finite bound

    p_t(x) ≤ binom(32,3)/[π⁵(8√6)³ Rδ⁶]                     (6)

is valid off the δ tube; no triple enumeration is required. It is intentionally a qualitative domination, not a competitive covering budget. The global Fourier peak is O(t^−6) at large time. Fixed-time continuity and dominated convergence now give finiteness and continuity of G on every compact set avoiding S2.

Consequently G attains a minimum m on D. Signed-permutation symmetry and −5 homogeneity give G(ry)=r^−5G(y)≥G(y) for 0<r≤1; G(0)=∞. Payment (4), the A038 upper obstruction, and compact continuity imply the exact identity

    essential inf_[-1,1]^6 G = min_D G = m.                    (7)

The essential upper inequality is legitimate even if a minimizer lies on a face edge or tie: continuity off S2 gives a relative cube neighborhood of positive six-dimensional measure with values arbitrarily close to m. This is a genuine finite-domain reduction with an actually paid global outer layer, but not yet a practical lower certificate on D.

## 4. A cutoff-safe continuous checker and its exact limitation

Let B_n(x)=∫_(1/n)^n p_t(x)dt, n=2,3,… . It is smooth everywhere, including S3, and increases strictly to G at every point of D. By Dini's theorem the convergence on D is uniform, using the continuity just proved. Therefore:

* If m>g*, one finite n exists with B_n>g* by a uniform positive margin on D.
* That fixed B_n has finite positive time moments and matching band gradients. A015's band Taylor certificate applies on every convex cell, including cells meeting S3 and face boundaries; no derivatives of full G are used.
* Compactness then gives a finite covering by passing neighborhoods. A concrete checker must retain an exact whole-root cover, and every leaf must be either a paid S2 tube leaf or a certified band leaf. A same-plane tube witness for all vertices is sufficient; choosing a different closest plane at every vertex is not. A fair rational subdivision whose diameters tend to zero is a conditional termination mechanism, not a claimed affordable run.

On the boundary d2=δ0 the larger paid value (4) leaves a strict margin; uniform band approximation therefore handles cells crossing that boundary when a same-plane tube test alone does not pass. This avoids silently leaving a curved seam uncovered. There is no volume tolerance for uncovered leaves.

This also exposes a precise endpoint obstruction. If m=g*, every finite B_n is strictly smaller than g* at a minimizer, since both omitted positive time intervals have strictly positive mass. No purely positive-band lower certificate of the form G≥B_n can then certify the exact floor at that point. It needs a strict global margin or an additional rigorous positive omitted-tail payment. If m<g*, the target is false. Current data do not distinguish these alternatives.

The existing full-G lower values do not supply B_n lower values: G≥B_n has the opposite direction for that inference. Likewise full gradG cannot substitute for gradB_n. Matching H/band values and gradients, or verified upper bounds for omitted masses, are new information required by the checker and were not computed in this round.

## 5. What the saved gradients actually imply

For each A038 point c, the exact interval dot product gradG(c)·(w−c) is negative for the same feasible ordered-face vertex w=(1,1,1,1,1,0). The fixed exact checks preserve all six components, with x1 displacement zero. At these supported generic differentiability points, each is therefore not a local minimum of the ordered face.

This supplies no uniform negative gradient bound on a neighborhood, no sign after reaching S3, and no region on which a flow is paid. Three pointwise derivatives cannot delete an open cell or prove a boundary-location theorem. In particular no conclusion min at first4=1 follows.

## Conclusion

New proved information is the S3 divergence obstruction (1)–(3), a precisely paid compact target domain (4)–(7), and a cutoff-safe conditional finite-certificate interface with its equality obstruction. There is no proved dimension drop or quantitatively affordable continuous cover. The global floor 1/65000 remains unresolved. No new kernel was used.
