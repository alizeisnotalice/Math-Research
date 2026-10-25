# A007: local J2 curvature certificate; radius improves, full-face gap remains

The finite batch completed normally in **205.325945708 seconds**, Python PID
**17359**, execution session **91632**. Eight new H_tau nodes and five real
gradient evaluations were computed at one fixed minimum-near center. No
round006 boxes were rerun, no MATLAB was called, and no commits or pushes
were made. The batch is stopped.

The best accepted clipped positive-face box has epsilon=2/125=0.016 and
radius **717/50000=0.01434**. Its continuous certificate is

    G(x) >= H_(2/125)(x) > 16666817/10^12 > 1/60000.

This is a local model result. It does not certify global g or any new C6
bound. The general L1/all-positive-scale operator interface and the integral
I certificate remain primary integration responsibilities.

## Inputs and actual arithmetic

The center is the previously exact-G-certified point

    (1000000007,998836847,892738545,807464540,467291939,96968636)
      /1000000007.

The old source is round003 `local_minimum/exact_local_reference_results.json`.
A read-only scan of main-clone upper JSON records found no existing H_tau
values at this exact center, so the eight requested tau nodes are new:
0.001,0.002,0.004,0.008,0.016,0.032,0.064,0.128. Existing G was not rerun.
The old center source, fan, F006 report and moment results, and local proof
inputs are SHA-256 locked in `input_manifest.json`.

Arithmetic was python-flint **Arb384**. Every H/gradient evaluation represents
the full **967680** signed-permutation/fan terms. Processing one member of
each global-sign pair with multiplicity two is exact. Gradients are real
automatic differentiation of the atan divided-difference recurrence, with
the copied-node DD6 identity for any repeated nodes. No complex-step method
was applied to a real-part formula.

Before the kernel, distinct-node gradients were independently checked against
both the copied-node identity and the differentiated barycentric formula;
a confluent anchor was also checked. Those actual scalar results are in
`anchor_calibration.json`. The five same-center gradients use tau<=0.016.
Their H values agree with the earlier eight time values; these shared AD
values are not new spatial nodes. Counting each gradient pass as another
evaluation gives 13 passes, below the 16-node budget, with one center only.

## Continuous moment and curvature certificate

For each band [a_k,a_(k+1)], use the nonnegative mass
Delta_k=H_(a_k)(c)-H_(a_(k+1))(c). Both bounds are rigorous Arb endpoints;
positivity permits intersecting with [0,infinity). In this batch every band
has a strictly positive lower endpoint.

The initial computation paid each band exp(r/a_k), not exp(r/epsilon). F007
then independently proved the stronger real comparison

    Phi(u)=exp(12 asinh(u/12)),
    p_t(c+h) <= Phi(||h||_infinity/t) p_t(c).

The proof uses the exact one-dimensional Cauchy ratio, a cube-facet
Caratheodory representation with at most six nonzero coefficients, and
Jensen's inequality. Its read-only source is F's `sparse_shift_bound.md`,
SHA-256 **c386ed307d9df8379c4dcca9983bf1571405b112b093c4bb9c9c0cd1dc04fede**.
`postprocess_sparse.py` applied it to the completed band data with **zero new
kernel nodes**.

Throughout the whole clipped box, the final bound is

    J2_epsilon(x) <= U2(r)
      = sum_(k=l)^6 a_k^-2 Phi(r/a_k) upper(Delta_k)
        + T^-2 Phi(r/T) upper(H_T(c)),  T=0.128.

The full tail [T,infinity) is included. The previously proved positive
curvature inequality becomes H_ii <= (4225/4224)U2(r). For a lower certificate,
coarea additionally gives the full-direction inequality

    D_h^2 H >= -||h||_2^2 J2/16.

The positive Taylor integral kernel then yields

    H(c+h) >= H(c)+grad H(c).h - ||h||_2^2 U2(r)/32.

The actual five clipped displacement intervals are used in the linear term,
including the short one-sided allowance toward x2=1; x1 is fixed. Every
center-to-box segment remains inside the convex clipped box. No
epsilon^-2 H replacement is used to accept these boxes. Details are in
`local_theorem.md` and the independently checked endpoint/rational fields
of `sparse_shift_results.json` and `checked_results.json`.

At the best box:

|Quantity|Rigorous Arb enclosure near the shown value|
|---|---|
|H_(.016)(c)|1.99863692360e-5|
|J2 at center, lower|0.00911922026054|
|J2 at center, upper|0.03650778367494|
|Whole-box J2 upper|0.08065620000403|
|Band [.016,.032] contribution|0.06603073911870|
|Band [.032,.064] contribution|0.01343416482697|
|Band [.064,.128] contribution|0.00115673004225|
|Tail [.128,infinity) contribution|0.00003456601612|
|Whole-box H_ii upper|0.08067529474835|
|Taylor curvature loss upper|2.59152907524e-6|
|Certified lower endpoint|1.66668174025324e-5|

The exact rational `16666817/10^12` is verified directly against the recomputed
Taylor enclosure. It avoids rebuilding a wide Arb ball from a short printed
radius: such radius strings have only a few significant digits and are not
the canonical acceptance fields. The checker reconstructs the formula from
the narrow H/gradient input balls and exports separate lower endpoints.

## Fair comparison and stopping decision

The accepted final radii are:

|epsilon|Time-band J2 Taylor radius|Matched uniform epsilon^-2 H Taylor radius|
|---|---|---|
|.001|.006560|.000840|
|.002|.007903|.001674|
|.004|.009435|.003325|
|.008|.011820|.006471|
|.016|.014340|.010758|

For the matched comparison, keep the same strict gradient, clipped geometry,
negative Hessian coefficient, and sparse spatial factor, and replace only
the J2 upper bound by epsilon^-2 Phi(r/epsilon)H_epsilon(c). This is a
diagnostic comparison, not a fallback in the accepted certificate. At the
best box that uniform budget is about 0.191147842894, whereas the retained
J2 budget is about 0.0806562000040, or **0.421957155168** of it.

The full method improves the best old single-time multiplicative radius
among all eight computed tau values from **0.002906** to **0.01434**, a factor
**7170/1453=4.93461803166**. Holding gradient/Taylor geometry fixed, the
time-moment contribution alone increases the optimized radius from 0.010758
to 0.01434, a factor **2390/1793=1.33296151701**. The sparse Phi improvement
over the original exponential-band Taylor certificate is only
**4780/4779=1.00020924880** at the optimum. The approximately fivefold gain
must therefore not be attributed entirely to time moments or to Phi.

This is below F006's suggested one-order-of-magnitude gate for expanding
coverage. At radius 10*0.002906=0.02906, the current sufficient Taylor
budget has upper endpoint strictly below 1/60000 for **all five tested
epsilons**; at epsilon=.016 it is about -2.47496443445e-6. This strictly
diagnoses failure of this allocated certificate budget at that radius. It
does not prove failure of the actual G bound. No grid enlargement is started.

## Coverage retained explicitly

All five accepted boxes have the same center and are nested; the best box
contains the others. Its five free-coordinate intervals are pairwise disjoint,
so its 120 positive-face permutation images are disjoint. Their exact total
volume is

    614593713570664440310497 / 488281253417968750000000000000
      = 1.25868791658190536e-6 approximately.

The exact uncovered volume is the complement of this fraction, recorded in
`checked_results.json`. The point `(1,19/20,9/10,17/20,4/5,3/4)` lies outside
all images. The complete positive-face complement remains unresolved.
Minus-five homogeneity extends covered directions inward in the cube only;
it does not resolve that complement.

## Reproduction and retained files

From this directory, use the existing Python environment:

    /Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python run_moments.py
    /Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python postprocess_sparse.py
    /Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python check_results.py

The first command is the bounded kernel batch; the latter two are pure
postprocessing. `postprocess_sparse.py` requires F007's exact-hash proof at
the recorded source path. Minimal artifacts are the registration, theorem,
this report, three Python scripts, input manifest, scalar anchor calibration,
raw results, sparse-shift results, and checked results. Terminal transcripts
and Python caches are unnecessary and removed after completion.

The local numerical evidence is real and the acceptance arithmetic is
replayable. The postchecker does not claim a second complete kernel replay;
F independently audits the mathematical interfaces and scalar reconstruction.
No global g/I or C6 advancement is entered into the ledger by this subtask.
