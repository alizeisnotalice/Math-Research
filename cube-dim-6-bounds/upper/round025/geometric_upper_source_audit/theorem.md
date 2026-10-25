# C025 geometric statements and their limits

Throughout, cubes are centered and axis parallel in R^6, Lebesgue measure is
used, and Q(x,r)=x+[-r,r]^6 has a common radius r in all coordinates.
All statements in Sections 1–3 below are proved here; Section 4 is explicitly
an unproved acceptance target. No global improvement of C6<=64 is obtained.

## 1. A uniform overlap improvement is false

For 0<e<1/2 and s in {+1,-1}^6 put

    c_s=(1-e)s,       Q_s=c_s+[-1,1]^6.

These 64 equal-radius cubes contain 0 in their interiors. Distinct centers
have sup-norm distance 2(1-e)>1, so no center lies in any other cube. In
particular they are a Besicovitch family in the center-exclusion sense.
Direct Cartesian products of one-dimensional intervals give

    intersection_s Q_s = [-e,e]^6,
    |Q_s intersection Q_t| = 64 e^k  if s,t differ in k coordinates,
    union_s Q_s = [-2+e,2-e]^6,
    sum_s |Q_s| = 4096.

Consequently the common intersection relative to one cube is e^6, every
pair overlap ratio is at most e, and

    (sum_s |Q_s| - |union_s Q_s|) / sum_s |Q_s|
       = 1-(1-e/2)^6 -> 0.

Thus neither a positive common intersection, pairwise center exclusion,
nor a fixed cardinality of 64 yields a uniform positive normalized volume
overlap. This refutes each proposed uniform-overlap refinement even for
equal radii and interior common intersection.

This does NOT show sharpness for the fixed-Lebesgue maximal operator:
centers alone form a zero-volume set, and no L1 input and positive-volume
bad set attaining ratio 64 have been constructed. General-measure
Besicovitch sharpness is not used to infer Lebesgue sharpness.

## 2. A genuine local six-dimensional strict gap

**Proposition (single support cube and bounded radius).** Let Q=Q(c,R),
R>0, let f>=0 be integrable and supported in Q up to a null set, and let
lambda>0. Assume M=integral f >= lambda |Q|. Define

    E={x: sup_{0<r<=R} (2r)^(-6) integral_{Q(x,r)} f > lambda}.

Then, with delta=43/15625 and gamma=64-delta,

    lambda |E| <= gamma M.

This proposition is an L1 statement for every positive R, with no
assumptions of continuity, boundedness, or symmetry on f. Its all-radius
consequence under the same support/height condition is proved below.

**Proof.** Translate and dilate to Q=[-1,1]^6 and R=1, retaining the
transformed height and mass, and put a=M/(64 lambda)>=1. All bad centers
lie in 2Q up to its boundary: outside it, a radius at most one cannot meet
Q in positive measure. Hence |E|<=4096.

Take eta=2/5, t=(1-eta)^6=729/15625. For each sign vector s, let

    C_s={y: 1-eta<=s_j y_j<=1 for every j},
    B_s={x: 2-eta<=s_j x_j<=2 for every j}.

The 64 input corner caps C_s are pairwise disjoint; the 64 center boxes
B_s are pairwise disjoint subsets of 2Q and each has volume eta^6.
If x in B_s has a witness with r<=1, then r>=1-eta, and
Q(x,r) intersection Q is contained in C_s. Therefore

    integral_{C_s} f > lambda [2(1-eta)]^6 = 64 lambda t.

The strict inequality follows from the strict bad-set threshold. The
case r<1-eta has empty intersection; the case r=1-eta can only strengthen
this necessary condition. If a<=22t, at most 21 caps can exceed
64 lambda t, since 22 such caps would have total mass strictly greater
than 64 lambda 22t>=M. The other 43 boxes B_s are disjoint from E. Thus

    |E| <= 4096-43 eta^6,
    lambda|E|/M <= (64-43 eta^6/64)/a <= gamma.

If a>22t, the crude support bound instead gives

    lambda|E|/M <= 64/a < 64/(22t)=500000/8019 < gamma.

Finally eta^6=64/15625, proving gamma=999957/15625. The changes of variables
multiply |E|, |Q|, and M by the same spatial Jacobian, so the original R
and center are restored. The cube integral is jointly continuous in (x,r) for r>0 by dominated
convergence: indicator functions converge away from cube boundaries,
which have Lebesgue measure zero, and the integrable dominator is f.
Consequently the maximal set is measurable. QED.

**Corollary (all radii at a support-compatible height).** Under the same
support hypothesis M>=lambda|Q|, the conclusion also holds for the full
centered cube maximal operator M_c over every r>0.

Indeed enlarge Q concentrically to Q'=Q(c,R'), where
R'=R (M/(lambda|Q|))^(1/6)>=R. Then f is still supported in Q' and
M=lambda|Q'|. For any center and r>R', the cube average is at most
M/(2r)^6<lambda. The full bad set therefore equals the one with r<=R',
and the Proposition applies. In particular this covers every compactly
supported input at all heights no larger than its mean over a containing
cube. It does not cover arbitrary heights for that input, nor all L1
inputs. This analytic extension uses no additional scalar checks.

The center boxes and input caps use all six coordinates simultaneously;
this is not a claim that the planar Lin bound tensorizes. The constants
were not optimized. The only evaluated geometric parameters were e=1/100
and eta=2/5 in eight exact rational sanity checks.

## 3. Where the covering proof needs new information

For f>=0, lambda>0, let selected cubes Q_i be disjoint up to null sets,
with integral_{Q_i} f >= lambda |Q_i|, and suppose a measurable bad set E
is covered by V=union_i (1+tau)Q_i, tau>1. Write

    M=integral f, S=sum_i |Q_i|, a=(1+tau)^6,
    B=sum_i (integral_{Q_i} f-lambda|Q_i|),
    O=integral_{(union_i Q_i)^c} f,
    D=a S-|E|.

All are finite when M<infinity, and B,O,D>=0. Disjointness gives
M=lambda S+B+O, and hence the exact identity

    a M-lambda|E| = a(B+O)+lambda D.                 (1)

The volume loss D includes both enlarged-cube overlap and points covered
by enlarged cubes that are not bad. The source proof estimates each term
on the right merely by zero. Section 1 forbids deducing a uniform gain
from the original-family overlap alone. A proof using (1) must exploit
its coupling with the input mass, not just cube geometry.

The all-radius Corollary applies to each g_i=f 1_{Q_i} at the original
threshold, but cannot be summed to control the full maximal bad set:
a witness may collect mass from several cubes or unselected input, and
can be bad while every component average stays below lambda. The
all-radius Corollary removes an isolated-component radius restriction;
the foreign and unselected mass estimates remain genuinely missing.

## 4. A precise next acceptance target (UNPROVED)

Put theta=delta/256=43/4000000. For each f>=0 in L1 and lambda>0, seek a
countable disjoint family Q_i=Q(c_i,R_i) satisfying
integral_{Q_i} f>=lambda|Q_i| and the following leakage estimate. Define

    g_i=f 1_{Q_i},
    L_i={x: M_c g_i(x)>(1-theta)lambda},
    F={Mf>lambda} minus union_i L_i.

The target is

    lambda|F| <= (delta/4) integral f = (43/62500) integral f.  (L6)

This is a concrete selection-and-mass lemma, not a consequence established
by the literature audit. It may be false. Even verification for finite
rational step-function configurations would only be a finite model;
passing such tests would not prove (L6) for general L1 inputs.

If (L6) were proved in the stated generality, the all-radius Corollary at
height (1-theta)lambda and disjointness would give

    lambda|{Mf>lambda}| <= [gamma/(1-theta)+delta/4] integral f
                          <= (64-delta/2) integral f.

The last algebra follows from
(64-3delta/4)(1-delta/256)=64-delta+3delta^2/1024.
Thus the explicit prospective constant would be 64-43/31250, but it is
**conditional and is not a new bound**. The next acceptance decision is
to prove or refute (L6), including foreign and unselected mass; there
is no authorization here to launch that work, a grid, or optimization.
