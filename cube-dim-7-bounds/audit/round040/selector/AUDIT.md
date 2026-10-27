# Independent audit: R040 fresh-witness-mass selector

**Verdict:** PASS for the claims actually made in the frozen proof. No blocking mathematical defect was found. The finite selector-to-halo mapping, its passage to arbitrary L1/all positive radii, the crude exponential estimate, and the stated many-source obstruction are valid. This audit does not supply the missing improved uniform halo bound.

**Frozen source:** `/Users/zhengzhihao/Desktop/T/cube-r040-cover-selector/cube-dim-7-bounds/upper/round040/selector/PROOF.md`.

**Source SHA-256 verified before the audit and by the independent microcheck:** `19d057f510e0f8b529f1090014cfc5a56f85efb9cd72fe6bdc28bfffe6b11b40`.

The proof was read in full, without reading or running its old verification script. The independent check below is a new parameter case, not a rerun of its three baselines. The source was not edited. The audit used the math-b01 finite certificate/normalization discipline, with direct self-contained reasoning; no external theorem beyond basic L1 integration/regularity was invoked. No MATLAB, production scan, fetch, rebase, commit, or push was used.

## 1. Selector order, frozen ownership, and original witnesses

The order is imposed on the **entire finite input family**, in nonincreasing witness radius with fixed ties. Therefore every previously selected source at a rejection has R_j >= R_i. This is enough for the later causal scale test; selected radii alone need not encode the rejected positions.

The selected F_j are the genuine differences `V_j minus union(previous selected V)` and partition the final selected union H, up to irrelevant faces. They are pairwise disjoint and the fresh selection inequality yields

    sum_j integral_(F_j) f <= ||f||_1,
    theta lambda sum_j |V_j| <= ||f||_1.

For the datum `g=f 1_H`, restriction does not change **any selected** witness integral: each selected V_j is contained in the final H. Also `g 1_(F_j)=f 1_(F_j)`. Thus both admissibility inequalities (4), including the strict original selected-witness test, hold for one common input g. This step does not replace the source masses by independent variables.

At a rejection the fresh inequality fails strictly (`fresh < theta lambda |V_i|`), while the original witness is strict. Subtraction gives the strict captured-mass lower bound (3). The selected prefix is exactly the union of the frozen regions created by then. Summing all final source pieces with R_j >= R_i only increases the mass compared with that prefix; future equal-radius sources do not invalidate the inequality.

Conversely, admissible source-only data satisfy the fresh selection test in their specified order because their F_j are precisely the previous-selected differences. No circular selector construction or altered witness is used.

## 2. Causal-spill definition and normalization

For any test radius R, `sum_(R_j>=R) g_j <= g` pointwise a.e. The upper cutoff is R_j, the radius of the **original witness**, and the averaging cube has radius `(1+1/a)R`. These are different roles and the frozen proof keeps them distinct.

Let c_a=1+1/a. For x in the small cube `Q_i=Q(y_i,R_i/a)`, every point of V_i is at distance strictly less than c_a R_i from x. Hence

    V_i subset Q(x,c_a R_i),
    |Q(x,c_a R_i)|=c_a^n |V_i|.

For a rejected item, its causal active sum contains the entire mass of the selected rejecting prefix inside V_i. Thus the test R=R_i has average strictly greater than `(1-theta)lambda/c_a^n`.

For a selected item, its own source piece is active at R=R_j and contributes at least `theta lambda |V_j|`. The same average is at least `theta lambda/c_a^n`, which is **strictly** above the halo threshold precisely because theta>1/2. Selection at equality in the fresh test causes no strictness problem.

Therefore selected and rejected small cubes are both covered by the **same** strict halo. The argument is pointwise for every x in each open small cube; no assertion that the original small-cube union was already an original maximal level set is needed.

The definition of Psi has the right numerator and denominator. Applying it to the actual datum gives `lambda |U| <= Psi ||g||_1 <= Psi ||f||_1`. Normalizing g by lambda normalizes all admissibility and halo inequalities simultaneously. There is no extra factor of lambda or c_a^n to restore. Translation and simultaneous spatial homothety of the datum and its source radii preserve this ratio.

The finite source cutoff changes only at finitely many R_j. At an endpoint R=R_j the active list uses `>=`; approximation from **below** retains those sources. Continuous dependence of moving-cube integrals on center and positive radius then justifies the stated countable-radius measurable version. Approximation only from above would not suffice, but the frozen proof explicitly allows the needed one-sided approximation.

## 3. Finite family to general L1 and all positive scales

The compact-level-set reduction is complete. Given any original strict level-set point x, choose a genuine positive witness radius R_x and the centered open small cube `Q(x,R_x/a)`. These cubes cover a compact subset of that level set and have a finite subcover. Each selected witness in this finite family still satisfies its original inequality for the same arbitrary input |f|. The finite selector theorem applies to that exact family.

The level set is open because for each fixed R the cube average of an L1 function is continuous under translation. Lebesgue inner regularity then passes the compact bound to the full level set, also if its measure was initially allowed to be infinite. No compact support, radius cutoff, rational scale lattice, or uniform bound on the number of chosen witnesses is required.

Thus `C_n <= Psi_(n,a)(theta)` holds for every fixed a>1. The liminf assertion follows because this lower inequality holds at every a, not from a limiting selector or a convergence claim for data.

## 4. Crude bound and circularity check

If T is positive at x, some active source contributes positive mass at a radius R <= R_j. Its support lies in V_j, so positivity of the intersection gives

    ||x-y_j||_infinity < R_j+c_a R <= (1+c_a)R_j.

Consequently the positive support of T is contained in the union of the original source cubes dilated by `1+c_a=2+1/a`. The union bound and the genuine fresh budgets give

    lambda |halo| <= (2+1/a)^n lambda sum_j |V_j|
                  <= (2+1/a)^n ||g||_1/theta.

This is a noncircular proof that Psi is finite, independent of the unknown value of C_n. Combining it with the interface and taking a to infinity, then theta to one through values below one, proves the claimed crude `C_n<=2^n`. No endpoint datum at theta=1 is needed.

The separate inequality `T<=M_n g` gives exactly `Psi<=c_a^n C_n/(1-theta)`. Its factor is larger than one, so combining it with the interface cannot improve C_n. The frozen proof identifies this correctly. No mass contraction follows from `g=f 1_H`: equality of the total masses can occur. The claimed improvement remains open, rather than hidden in this recurrence.

## 5. All-dimensional many-source counterexample

The source claims this obstruction for every n>=2. That scope is correct; no n=1 extension is proved or needed for the stated claim.

The V_j are genuine centered cubes: their centers are `(delta_j-R_j,0,...,0)` and every coordinate half-width is R_j. The order is strict because `R_j=3-j/(2N)` decreases with j, with minimum 5/2>1, so every source precedes the target radius-one cube.

For every point in P_j or B_j, its first coordinate obeys

    delta_(j-1) < x_1 < delta_j.

Every earlier source has right face at most delta_(j-1), so the point is **strictly outside** every earlier source. The current left face is negative (`delta_j-2R_j <= 1/2-5 < 0`), while x_1>0 and x_1<delta_j. Its other coordinates have absolute value less than one, whereas R_j>=5/2. Thus each new support rectangle lies **strictly inside the current open cube**. Endpoint coincidences between open rectangles and cube faces do not create positive-volume boundary exceptions; the pointwise inequalities above are strict. No compact containment of their closures is required.

An earlier support rectangle lies entirely inside each later source; a later support rectangle lies entirely outside an earlier source. The target support P_* has first coordinate >1, whereas all source right faces are <1/2, so it is outside every source. The positive-mass support rectangles are pairwise disjoint.

Every bulk mass is strictly positive, uniformly in the allowed parameters, because

    2|V_j|-A/N >= 2*5^n-(1-theta/2)*2^n/N > 0.

The actual fresh region of source j therefore carries exactly its own two rectangles, of total mass 2|V_j|. Its full witness contains all earlier rectangles and its own, hence has mass `sum_(ell<=j) 2|V_ell| > |V_j|`. This verifies strict original witnesses and all source decisions for one common f.

Inside the target, every P_j is fully present and every B_j is absent because its second coordinate is negative; this is where n>=2 is used. P_* contributes its full B outside the selected union. Directly,

    (A+B)/v_* = 1+theta/4 > 1,
    B/v_* = 3theta/4 < theta,
    A/v_* = 1-theta/2 > 1-theta.

Thus the final witness is strict and genuinely rejected, and the frozen source pieces each contribute exactly A/N. For any proposed fraction alpha>0 of the rejection threshold, choosing

    N > (1-theta/2)/((1-theta)alpha)

strictly defeats that claim. The same construction defeats a fixed fraction of the total captured mass by the exact share 1/N. No input-mass normalization or difficult operator ratio is claimed; the large bulk input is retained and counted.

## 6. Independent new exact microcheck

`microcheck.py` was actually run with Python 3.15.0 standard-library Fraction arithmetic on

    n=4, N=3, theta=5/8, lambda=7/5, a=9/4.

This differs from all three frozen baseline cases and also checks a nonunit level and another a. It independently computes support intersections, true witness masses, fresh regions/decisions, actual target contributions, halo-threshold normalization, and total mass conservation. No original baseline script was run.

All assertions passed. The normalized target values were:

| Quantity divided by lambda times target volume | Exact value |
|---|---:|
| Total witness mass | 37/32 |
| Fresh mass | 15/32 |
| Captured mass | 11/16 |

Each source contributes exactly one third of the captured mass and `11/18` of the rejection threshold. These numbers check this one finite case only; the all-n>=2 conclusion rests on the inequalities in Section 5, not on a numerical extrapolation. `microcheck.json` records the full exact masses and decisions. Reproduction: run `python3 microcheck.py` in this audit directory.

## 7. Integration limits

No correction to the frozen proof is requested. An integrating report should retain these boundaries:

* The uniform estimate for Psi that improves n log n is unproved.
* The single-source counterexample is all n>=2, not n=1, and is a selector-allocation obstruction rather than an operator lower bound.
* The common-spike discussion assumes each witness contains the entire input support. It does not establish a halo bound for all common-intersection families or all common witnesses.
* The source definition uses support in H in the a.e. L1 sense; cube faces and values on null sets do not affect the proof.

Only this independent clone's `cube-dim-7-bounds/audit/round040/selector/` directory was created for the audit. The parent retains proof and integration responsibility.
