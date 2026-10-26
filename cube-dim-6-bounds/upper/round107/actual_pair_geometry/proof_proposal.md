# 107 prepared: k1 exact pair geometry and a coefficient-free mixed-q budget

Status: proof/source proposal only. No107 pair generation, numerical pilot,
cell evaluation, kernel, peak or old tuple/grid program has run. Two finite
configurations require a separately hash-bound parent release. The user
runtime policy effective2026-10-26T16:43:45.633722UTC removes previous
program wall/CPU caps; all old actual receipts remain unchanged. This
proposal instead has fixed inputs,1024 leaves/2040 cell bounds/depth20,
checkpoints and explicit mathematical stopping conditions.

The target remains dimension6, fixed Lebesgue measure, centered axis-
parallel common-radius cubes, general L1 inputs and every positive radius.
The original equal-weight body32 alpha1/common-time Cauchy frame and
Fourier normalization are unchanged. Only k1 is included. No lower-bound,
Astra, MATLAB, agent expansion, stage or push is proposed.

## 1. Actual pair coefficients, with no geometric-orbit assertion

Let A={1},B={2,...,6},v=ones,w=(1,-1,-1,-1,-1,-1). Body representatives
are exactly(1,s2,...,s6), each sign tuple once. For each retained b=v or w,
visit every ordered pair(i,j) with i!=j and i,j!=b:930 visits per branch.
There is no additional unoriented sign factor. Face coordinate1 stays+1.

Write Bi=sum_(ell=2)^6 i_ell, Bj similarly, c_ij=i.j and Delta=36-c_ij².
For x0(u)=(1,u,u,u,u,u), the true pair-plane center distance is

 P_ij(u)=1+5u²
  -[6(1+Bi u)²-2c_ij(1+Bi u)(1+Bj u)+6(1+Bj u)²]/Delta. (1)

The coarea is J_ij=sqrt(Delta), with Delta in{20,32,36}. The first
configuration aggregates only equal exact triples of rational quadratic
coefficients and equal Delta, retaining every ordered multiplicity n_C.
This is algebraic(P,J2) aggregation, not a claim that its classes are
geometric orbits. Each branch must retain930 visits and determinant
counts180,450,300. Bi,Bj have6 values and c_ij5 values; after their swap
symmetry there are at most21*5=105 algebraic classes per branch.

Every quadratic minimum/maximum on[1/4,3/4] is certified exactly from its
two endpoints and, when present, the rational vertex. Each minimum must
be strictly greater than delta², delta=1/64. Thus each individual

       sqrt(P_C(u))-r>0, 1/4<=u<=3/4,0<=r<=delta,         (2)

has its own positive certificate. The051 all-other-plane theorem already
provides a stronger existence proof: a pair avoiding b cannot be the
plane span(v,w), whose only body lines are v,w. It therefore stays away
from this tube. That uniform theorem is used for admissibility only;
the cost never substitutes its global g_min for the individual P_C.

## 2. Real clamp functions for the restricted source events

For k1 the true face coordinates are x=(1,u1_B+z), sum z=0, r=|z|.
The perturbation from x0(u) is exactly r. The first-two-distance source
gaps may be chosen as

 hv0=sqrt(min_(|a|<=1/2)[(1-a)²+5(u-a)²])-r,
 hv1=sqrt(5/6)min(u,1-u)-r,
 hw0=sqrt(min_(|a|<=1/16)[(1-a)²+5(u+a)²])-r,
 hw1=sqrt(5/6)(7/8-u)-r.                               (3)

The origin minimizers are clamp((1+5u)/6,-1/2,1/2) and
clamp((1-5u)/6,-1/16,1/16). The corresponding squared-origin gaps at
r=0 have the exact pieces

 Fv=(5/6)(1-u)² for u<=2/5,
 Fv=1/4+5(u-1/2)² for u>=2/5;
 Fw=(5/6)(1+u)² for u<=11/40,
 Fw=(17/16)²+5(u-1/16)² for u>=11/40.                   (4)

The line minimum switches at u=1/2. Hence all initial u intervals have
endpoints1/4,11/40,2/5,1/2,3/4. Further bisections stay inside one
valid regime. The existing051 weaker gaps prove positivity everywhere;
the cell code also verifies every square-root/linear lower gap before
division. No(1-20r)^-5 or other global relative-gap substitute is used.

For each retained branch define

 W_b(u,r)=sum_C n_C/[sqrt(Delta_C)(sqrt(P_C(u))-r)],
 H_b(u,r)=h_b0(u,r)² h_b1(u,r)²,
 S(u,r)=W_v/H_v+W_w/H_w.                               (5)

All930 pairs remain in each branch. Every source event and normal
cancellation is retained through104's actual restored four-dimensional
quotient; enlarged small events are used only as positive upper bounds,
with no independence or disjointness assumption. Only two deletions and
restoration precede the time integral. No third-direction sum is added.

## 3. Factor the common integrated-density coefficient out completely

The accepted098 full-symbol coefficient gives

 Kold=512/(pi^4 c³)=3sqrt3/(64pi^4), c=32/sqrt3.

The104 time-first source interface yields

       G(x)<=Kold[1/(J_vw r)+S(u,r)], J_vw=2sqrt5.       (6)

The cells integrate only the coefficient-free shape. Thus a later
accepted common coefficient can be applied by pure algebra without
regenerating any geometry or cells. No pending105 or108 theorem is
silently used. A separately registered reference expression is

       Kref=1/(8pi^5), Kref/Kold=8/(3sqrt3 pi).           (7)

Its proposed origin is the parent's pending108 projection bound, not an
accepted input here. Every Kref number is explicitly a conditional shape-
expression diagnostic. The sole currently valid kernel upper certificate
uses Kold. Neither the sharp-Khinchine105 nor the108 claim is adopted as
a theorem or used to change the registered coefficients after a result.

## 4. Exact k1 spatial measure and removal of the origin singularity

For u in[1/4,3/4] and r<=1/64, all five coordinates u+z_i stay strictly
inside[0,1]. The face tube is therefore exactly the full product of this
u interval and a four-dimensional zero-sum ball. Its Jacobian is sqrt5;
there is no nontrivial maximal-coordinate cone in k1. If v4=pi²/2, its
radial measure is4v4 sqrt5 r³ dr du.

Define the fixed conservative shape integral

 Jshape=4v4 sqrt5 int_(1/4)^(3/4) int_0^delta
                   r³[1/(J_vw r)+S(u,r)]^(6/5) dr du.   (8)

This retains both enlarged complementary branches and the large branch
inside the q power, rather than splitting their mixed terms. Set
r=delta*s^5. Since3-q=9/5, exactly

 Jshape=20v4 sqrt5 delta^(14/5) int du int_0^1
         s^13[1/J_vw+delta*s^5*S(u,delta*s^5)]^(6/5) ds. (9)

The integrand in brackets is positive and continuous at s=0. No infinite
sample or unverified near-zero truncation is needed. The original center
line has zero face measure; the transformed finite integral pays its
whole neighborhood. For any common coefficient K,

       I_expression(K)=K^(6/5)Jshape.                   (10)

Equation(6) proves the actual k1 integral is at most I_expression(Kold).
A lower bound on Jshape concerns this fixed upper-majorant expression,
not a lower bound on true G, true I, the model's best cost, or C6.

## 5. Cell bounds are analytic ranges, not point samples

For each rectangle[a,b]x[s0,s1], put r0=delta*s0^5 and r1=delta*s1^5.
Compute each exact quadratic range P_C([a,b]) by its endpoints and
possible rational vertex. The monotone denominators in(5) then give

 n/[sqrtDelta(sqrt(Pmax)-r0)] <= term
               <=n/[sqrtDelta(sqrt(Pmin)-r1)].          (11)

Apply the same quadratic-range procedure to the actual squared-origin
piece(4), and exact affine ranges to the line gaps. Outward Arb384
square roots and strictly positive endpoint denominators enclose both
branch sums; all mixed terms remain inside S. This produces rigorous
S_lower<=S(u,r)<=S_upper on the whole rectangle. Then

 T_lower=(1/J_vw)_lower+r0*S_lower,
 T_upper=(1/J_vw)_upper+r1*S_upper                       (12)

enclose the regular bracket in(9). Fractional powers use outward Arb,
with additional exact rational endpoint checks power_lower^5<=T_lower^6
and power_upper^5>=T_upper^6. The monomial weight is integrated exactly:

 cell_weight=(b-a)(s1^14-s0^14),
 mass_constant=(10/7)v4 sqrt5 delta^(14/5).              (13)

Multiplying those positive endpoint bounds gives a valid lower/upper
enclosure of each cell's contribution to Jshape. The final exact rational
sum of leaf endpoints encloses(8)-(9). Refinement records a complete
parent/children tree from the eight initial rectangles; a volume sum
alone is not treated as proof of coverage.

## 6. One finite adaptive enclosure, two frozen comparisons

The initial cells are the four clamp-valid u regimes times the s intervals
[0,1/2] and[1/2,1], eight leaves. At every step refine the eligible cell
with largest absolute contribution upper-minus-lower. Split alternating
normalized u/s depth, so the choice is deterministic and does not add
sample probes. Exact weights determine the priority. The limits are1024
leaves,2040 total evaluated cell bounds and total depth20 from each seed.
Inputs, coefficient comparisons and caps are never refreshed.

Apply the same Jshape interval once to Kold and Kref. Let
T=64(g0/6)^(6/5)-I_paid_CEH, g0=1/65000 still unproved globally. If both
constructed expression lower endpoints exceed T_upper, stop with a
strict dual-expression obstruction. Otherwise continue to the fixed
node/depth cap and report the actual separation status of each expression.
The Kref obstruction/success remains conditional on its unaccepted
coefficient. This avoids stopping merely because the old coefficient is
too expensive when the independently proposed reference may still fit.

Compare the valid Kold upper endpoint directly with the saved104 k1
complete upper allocation lower endpoint, about1.795895285433e-4. A strict
inequality certifies a better upper-budget expression for the same full
k1 middle region. Passing T would pay only this k1 slot; every other
middle type, far and global g_L remain unpaid. A positive diagnostic
remaining allowance is not global S or a C6 improvement.

If both fixed expressions are obstructed, merely refining this numerical
enclosure cannot make either fit that allocation. This excludes only the
specified radial pair-distance/clamp/constant-coefficient expression. It
does not exclude actual angular distances, directional marginal bounds,
the weighted line-trace retaining two time denominators, a different
spatial majorant, or the true model. No automatic k2..5 extension follows.

## 7. Execution and local compute plan

The script is currently unrun and has no parent release file. Before
mathematical imports it requires exact script/registration release hashes.
It refuses existing runtime artifacts, checkpoints after configuration1
and every64 refinements, and stops on first failure without retry.

The latest user explicitly removed program runtime limits. This source
therefore contains no177-second fuse or CPU60 cap. It remains a finite
single-process first-feedback stage, with at most1860 new pair visits and
428400 aggregated-term cell evaluations by the105-class bound. Actual
cost is unknown before this first release; no timing pilot is run. The
main process can use other available cores for independent upper work.
Parallel cell execution, if later warranted by actual cost, needs a new
explicit finite registration; this run will not expand itself to fill CPU.

No107 numerical result, actual strict improvement, bound obstruction or
new operator constant is claimed before the authorized run.
