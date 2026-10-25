# A015: valid band-Taylor cell certificates, eight old-data calibrations

The positive-band simplex/cell lower certificate is valid, including the
factor1/32 for body32, a correct same-plane singular-tube rule, and exact
midpoint-tree coverage. The finite old-data calibration completed with
**four PASS and four insufficient lower bounds**. No new kernel, gradient,
MC, grid or MATLAB computation was made.

Fixed body32 T=.032 is already refuted by B010. During this task the main
agent also supplied the B015 obstruction at fixed T=.048, whose checked
record was read and hashed. Neither time remains a viable full-domain
target here. The valid inequalities and local calibrations do not undo
those counterexamples. No global m_T or new C6 bound is claimed.

## Correct lower certificate

Let B_(a,T)=H_a-H_T=integral_a^T p_t<=K_T, with0<a=t0<...<tn=T.
At a single verified center c, form nonnegative time masses
Delta_i=H_ti(c)-H_ti+1(c). For a convex cell with vertices v and center c
inside it, put r=max||v-c||inf and D=max||v-c||2^2. The whole-cell moment
upper bound is

    U2=sum_i ti^-2 Phi(r/ti) upper(Delta_i),
    Phi(u)=exp(12 asinh(u/12)).

There is **no infinite-time tail**. All mass coefficients are positive,
and the upper subtraction endpoints are used. Coarea with R=A^T/32 gives
the actual lower Hessian

    D_h^2 B >=-||h||2^2 J2_(a,T)/16.

Taylor's integral remainder pays a further factor1/2. The resulting exact
cell lower certificate is

    K_T >= L_B+min_vertex lower(gradB(c).(v-c))-D*upper(U2)/32,
    gradB=gradH_a-gradH_T.

No Hessian upper estimate is used to infer a lower estimate, and no
convexity of B, K_T or the ordered chamber restriction is assumed.
Every gradient/value belongs to the same frame, center and time endpoints.
Segments stay inside the convex cell, and squared displacement and linear
extrema are bounded at vertices. The self-contained proof is `theorem.md`.

For body32 augmented by beta axes, the same proof gives the Taylor loss
D*U2/(32+beta^2), because the Gram matrix is(32+beta^2)I. Phi remains valid
using only the original sparse body translation, with extra axes unchanged.
This is a theoretical extension only: no old body H/gradient or A009 barrier
is transferred to an augmented kernel.

## The actual eight calibrations

The sole center is the frozen A007 rational vector

    (1000000007,998836847,892738545,807464540,467291939,96968636)
      /1000000007.

Existing A007 supplies H and gradients at a=.001,.002,.004,.008, together
with the intermediate time values. Existing B010 supplies H_.032 and its
gradient at exactly that same center; its H overlaps the archived A007
H_.032 value. The band data use only old endpoints between a and T=.032.

The geometric calibration cell is the max-norm r box about c, clipped to
x1=1 and[0,1]^5. Its32 vertices remain in the ordered positive face, and
an exact rational convex-combination witness verifies c is inside. The
clipped x2 upper endpoint1 is used in the linear minimum. For these cells,
the exact maximum Euclidean squared displacement is5r^2.

The target for comparison with the old Phi radius is
m0=T/4000=1/125000=8e-6. All shown values below are lower-certificate
formula values, not estimates of the cell's actual minimum:

|a|r=.004 lower|r=.008 lower|
|---|---:|---:|
|.001|8.70143157169e-6 PASS|-1.13496088053e-5 insufficient|
|.002|8.80590055727e-6 PASS|3.63629634628e-6 insufficient|
|.004|8.88274132543e-6 PASS|6.32049223733e-6 insufficient|
|.008|8.71277463254e-6 PASS|7.44759742192e-6 insufficient|

The same four/four flags hold at the weaker bare thresholdT/4096=1/128000.
An insufficient lower formula is not a counterexample to the actual K value.
All four passing a choices certify the **same** r=.004 cell, not four
additional covering regions. The best tested bound uses a=.004.

The saved B009 pure-Phi radius **991/500000=.001982** was for the matched
target1/125000. The new passing tested radius1/250=.004 is larger by exactly
**2000/991=2.01816347...**. The old bare64 radius was a different value,
.002277; it must not be mixed with the m0 comparison. This batch does not
bisect a new optimal radius or add any unregistered a/r.

All acceptance comparisons use exact rationals: input Arb lower/upper
endpoints are converted to exact dyadics, vertex lower endpoints are
minimized rationally, U2 is rounded up, and the final Taylor expression
is a rational number. `results.json` contains the complete rational cell
certificate and each positive-band contribution. Wide serialized Arb balls
or a status flag alone are not used to accept a cell.

## Complete future cover acceptance contract

`checker_schema.md` specifies all required source/geometry/model fields,
exact acceptance checks and termination statuses. It is a contract for the
future global checker, not a claim that a global cover or kernel-proof reader
was executed in this turn. Its required conditions include:

* A verified source catalog with exact center, frame, time and gradient
  identities; a file hash alone does not prove a kernel value.
* The closed ordered-face root with six prefix-one vertices, exact longest-
  edge midpoint children, and a finite tree with no missing child or leaf.
* A TUBE leaf must have one fixed A012 canonical plane containing every
  vertex within the proved delta. Different closest planes per vertex are
  not accepted.
* A BAND leaf must pass the rational Taylor lower certificate. Every
  terminal leaf must be paid at a common m_star>T/4096. Tiny uncovered
  volume, unresolved signs or incomplete nodes never produce GLOBAL_PASS.
* An approved model/time obstruction returns EXCLUDED_MODEL_TIME; in
  particular the fixed body32 T=.032 and T=.048 are excluded.

Exact signed-permutation symmetry and
K_T(r y)=r^-5 K_(T/r)(y)>=K_T(y),0<r<=1, transfer a genuinely complete
ordered-face certificate to the full cube. This fixed-T argument is
separate from a false fixed-H_a homogeneity claim.

The algebraic T48 tube component remains true:
rho<=1/100000 gives K_.048>=K_.032>=3/250000>3/256000. It only pays that
neighborhood. B015's compact-region point obstruction prevents it from
supporting a full-domain T48 claim. A new time/frame would need its own
model-matched rigorous values and complete paid tree; none was proposed
or tested by this task.

## Two exact simplex sanity tests

The first test proves why per-vertex closest-plane coverage is unsafe:
all six root vertices lie in S2, while the root barycenter has exact
distance squared1/9. At delta=1/100000 no one canonical plane covers all
vertices. Thus the root is rejected by the actual same-plane rule.

The second test checks the root's longest edge is(0,5), its midpoint gives
two exact children, and the barycentric replacement identities prove their
union is exactly the parent. Exact volumes are1/120 for the root and1/240
for each child. The coverage proof is the barycentric identity, not the
volume sum by itself. No function value was evaluated in either toy test.

## Execution and complete remaining gap

`calibrate.py` completed with exit0, Python PID **32194**, in
**0.024190958 seconds** of script time (tool chunk07516e). It ran exactly
eight existing-data cases and two geometry tests. Source hashes were
registered before computation. After the later B015 steering, only the
context/exclusion documentation and metadata changed; the numerical cases
and coefficients did not change. Executed and final text hashes are retained.

There are zero new kernel or gradient calls, MC samples, MATLAB calls,
agents, commits or pushes. All old directories remain unchanged. The full
ordered-face root is still unresolved; the two overlapping local calibration
cells do not form a global partition, and no global cover was attempted.

The inequalities are reusable proof tools, not evidence that a refuted
global hypothesis survives. The finite batch is complete and stopped.
Reproduce the eight numeric cases and two exact toys with:

    /Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python calibrate.py
