# A076 report: variable-amplitude rank3 layers and finite-frame semiconcavity

The round supplies a complete analytic proof submission for local upper
semiconcavity of the Green function of every finite positive spanning Cauchy
frame in R6, outside its **true S2**. It includes the original body32 as a
formal specialization. The parent retains final proof acceptance. No actual
body32 semiconcavity constant, scalar vertex certificate, continuous full-face
I, global g cover, or Cn improvement has been produced.

The theorem text is unchanged from the version sent for parent/G076 review:

    theorem.md SHA256
    5b9bc17105fbf7b22be55049ab488d8aa7c4f47cdd10b2bc8ec18d7baf63b7b9

The parent reported that its full read and G076's preliminary full read found
no fatal gap. This is preliminary review feedback, not an invented final
acceptance. The mathematical proof and the finite checks remain separate.

## What is new relative to sealed072

A072's convex three-dimensional deficit norm and finite-overlap product
lemma are accepted read-only inputs. Its axis6 proof, axis fixture, overlap
geometry fixture, and abstract leading-cusp counterexample were not repeated.

The new proof obtains an exact local decomposition of the **actual** finite
frame law. A fixed smooth cutoff partitions scalar parameters by the rank of
their large columns. Rank-at-most-two terms vanish on the qualified ball;
rank-at-least-four terms retain four tail factors under all cutoff derivatives
and have a C2 integrated remainder. Each incident rank3 span U is grouped
using every original column in U. Its full in-plane density is
t³ a_U(y)+t⁵ r_U(t,y), with a_U>=0 and uniformly controlled derivatives. This
follows from the full three-dimensional fan, odd time parity, and the three-
independent-large-parameter density bound.

The remaining coupled-variable argument has three essential steps:

1. Remove outside cutoffs from the **unexpanded bounded amplitude**. A basis
   containing its one large quotient column gives an integrable C2 error.
   There is no weighted column in this step, so parallel quotient columns
   create no false two-column basis assumption.
2. For each nonzero quotient column b_j choose a basis containing j. Positive
   coarea gives weighted density bounds M_(j,k)<=C_(j,k)t^(k-3), k=0,1,2.
   A different fixed basis is used for differentiation. The quadratic Taylor
   remainder then has second-derivative integrand
   t³[t^-2 M2+t^-1 M1+M0]=O(1), uniformly on the full local ball.
3. The linear signed moment is treated with its necessary all-time linear
   subtraction. Its exact renormalized kernel is

       -3/(4pi²) integral_S² sign(b_j.theta)
                  (theta.z)|theta.z| / psi(theta)^5 dtheta,

   which is C1,1. The finite-time layer also includes the explicit linear
   term and a C2 upper tail. The proof does not discard either term.

This yields G=R-sum_U a_U N_U, R locally C1,1. Different rank3 spaces may
overlap away from true S2; no geometric-disjointness hypothesis is used.
The unique large-column span is the parameter-pattern label. A072's product
lemma then gives a finite distributional Hessian upper bound. The constants
are finite formulas once the geometric, basis, fan, and amplitude-bound data
are supplied; none has been evaluated for body32 here.

The proof explicitly retains all moving-fiber effects through fixed right
inverses, all cutoff derivatives, full fan normalization, and the signed
moment's small- and large-time convergence. It does not differentiate an
event indicator or use040's positive J2 divergence as an actual signed
Hessian conclusion.

## Whole-box meaning and remaining numerical work

For a qualified orthonormal box with half widths h_i, a certified constant
C and true scalar vertex lower bounds l_v give

    G(x)>=min_v l_v-(C/2)sum_i h_i²

through concavity of G-C|x|²/2. The sharper multilinear variance expression
is in theorem equation(25). An orthonormal five-dimensional face uses its
32 vertices. A nonorthonormal parameterization must pull back the quadratic
form. The theorem lists all required geometric, quotient, complete in-plane
fan, remainder, and model-identity data for a future certificate.

The main remaining production work is to construct useful certified constants
for actual body32 neighborhoods and match them to actual scalar vertices and
the intended continuous-I contract. No affordable constant, entire face cover,
S2-crossing box, or global g replacement is implied by the qualitative proof.
No new assertion that body32 has divergent Hessian on S3 is made.

## Registered checks and actual results

Exactly three new small configurations were registered before execution,
under the task ceiling of eight. They all passed on the separately authorized
continuation, using local Python3.15.0 / SymPy1.14.0.

| Configuration | Actual work | Result and scope |
|---|---|---|
| C001 | Exact fixed-cutoff polynomial algebra | C2/C4 endpoint matching; derivative bounds2 and6; large-factor bounds1,4,20. No density evaluation. |
| C002 | One complete three-dimensional four-column density path, four upper-half-plane residues, six pair-span determinants | Frame [e1,e2,e3,(1,1,1)] at (1,3,7) is off its true S2. Exact complete residue sum is real/even after dividing by t³; positive leading coefficient and t⁵ correction are verified. |
| C003 | Exact weighted/signed-layer identities | Scalar weighted-sup identities, fixed-basis remainder exponents, positive sine Fourier phase, integration-by-parts identity, signed angular coefficient, and C1 matching of s|s|. No numerical moment integration. |

The new non-axis three-dimensional density has the exact result

    pi³ p_t(1,3,7)/t³
       =5(32t⁶+276t⁴+703t²+279)
        /[8(t²+1)(t²+4)(t²+9)(4t²+1)(4t²+9)(4t²+49)].

Its leading coefficient is155/14112>0; after subtracting this coefficient
the quotient by t² has finite value -920735/24893568 at zero. The six exact
pair qualification determinants are7,-3,4,1,-6,2. This is a tiny full 3D
symbolic density calculation, with four actual residues. It is not labeled
zero mathematical work, a 6D full kernel, or a real body32 certificate.

All three checks support limited algebraic parts of the written argument;
they cannot establish its uniform analytic claims by themselves. No old069,
old072, old052 axis, or old048 numerical reference was rerun.

## Startup failure, explicit continuation, and resource ledger

The original frozen guard set its own CPU soft/hard limit to2 seconds.
Its child inherited the hard limit and failed while trying to set27 seconds,
before importing SymPy or entering C001. The error was
`ValueError: not allowed to raise maximum limit`. Original guard PID85973,
child PID85974, tool chunk44b27b; guard elapsed wall0.037732709 seconds. This
invocation executed **zero** scientific configurations, residues, or pair
determinants. It was stopped and reported; its original source, registration,
and failed receipt remain unchanged.

The parent then explicitly released one new continuation with a separate
guard and supplement registration. The corrected guard uses soft2/hard29
seconds so that the unchanged checker can lower its own soft/hard limit
to27. It preserves the same three inputs and the first-scientific-failure
stop rule. It is not an unrecorded automatic retry.

Continuation guard PID86350, child PID86351, tool chunk772b26, both completed
normally with no live child left. Guard elapsed wall0.429159042 seconds;
the checker measured wall0.368355500 and CPU0.282006 seconds. OS child
user+system CPU was0.337961 seconds, including startup outside the checker's
timer. Across both invocations, recorded guard CPU plus child OS CPU totals
0.456506 seconds. Both invocations were under the180-second whole-program
limit. The guard's normal timeout path allows165 seconds for work and10
seconds to kill and wait its own process group, reserving further overhead
inside the180-second hard alarm. No timeout occurred.

Final counts: one zero-science startup failure plus one explicitly authorized
scientific execution; three new configurations, one new3D density path,
four residues, six pair determinants; zero body32 full kernels, zero6D full
kernels, zero old baselines, zero MATLAB calls, zero recursive agents, zero
automatic retries, zero staging/commits/pushes. No additional execution is
authorized by this seal. The guards refuse reuse of existing receipts/results.

## Artifact and provenance boundary

Only this076 directory was written. Accepted072 theorem and manifest were
read and SHA256-checked against the parent's supplied pins. Original069 and
other frozen work remain read-only. `registration.json` pins the original
checker, guard and input manifest. `supplement_registration.json` pins the
new guard and the original failed receipt, binding the explicit continuation.
`results.json`, both execution receipts, and the tool receipt preserve actual
runtime and process provenance. `manifest.json` seals the complete artifact
set; it is not a final mathematical approval by the parent.
