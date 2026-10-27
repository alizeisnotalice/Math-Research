# Final manuscript mathematical review

**Final status: ACCEPT; all three minor definition notes are resolved.**
The final source SHA256 is
`94c5e61d3bf86b16736e328ffe3e412d46df957dd2b8f8760939edc8831e8692`.
Its exact incremental check is recorded in Section 10. The full first-pass
review and its original line numbers are retained below for traceability.

2026-10-28. **ACCEPT for mathematical content, with three minor definition
notes below.** No substantive error was found in the new operator theorem,
its explicit threshold, the actual-history transfer, the older comparison
bounds or the stated pointwise-load obstruction. No numerical threshold,
constant or main proof needs changing. The three notes would improve the
standalone final presentation; none affects an application made in the paper.

## 1. Reviewed version and boundaries of this review

The complete 1542-line TeX source was read, with line numbering from
`/Users/zhengzhihao/Desktop/T/cube-cover-manuscript/cube-dim-7-bounds/paper/main.tex`.
Its reviewed SHA256 is
`e52394e0a80eb01338e38a2cddc7c9b25357f141433a2214b8cd25e7f7b6db37`.
All line references below refer only to that exact version.

The source was compared with the locked R042 full proof and explicit
threshold, and with the already independently reviewed R030-R037 sources.
Actual source identities are recorded in `sourcehash.json`. This review
does not merely accept a statement because a private source was audited:
the mathematical argument as written in the manuscript was checked from
its finite selection through the final scalar inequality.

This is a mathematical TeX review, not PDF layout or compiler validation.
No compilation, production run, old numerical rerun, n=7 calculation,
literature search, MATLAB or Git operation was performed. No manuscript
file was edited. Only this audit directory was written.

## 2. First-pass definition notes, now resolved

1. **Lines 407-412: define I_n(t) before first use.** Its explicit
   definition is currently delayed until lines 445-450. Insert, before
   the integral-uniform lemma, the sentence
   `记 I_n(t)=\int_1^t\ell_t(\rho)^{n-1}\,\mathrm d\rho/\rho` with
   `\ell_t(\rho)=1+1/(2t)-1/(2\rho)` (or move the existing definition
   there). The integral estimates themselves are correct.

2. **Lines 699-710: define the full source-shell measure explicitly.**
   The notation `\lambda^0_{\rm source}` first appears at line 708
   without an explicit density. Add, before the lemma, the definition

   ```
   d lambda_source^0(z)=2^(-n)||z||_infty^(-n)
     1_{t(1-3d/4)<||z||_infty<t(1-d/4)} dz.
   ```

   It would also be useful to say that each maximal face at fixed rho
   has radial weight d rho/(2rho), with n-1 free uniform coordinates.
   This is already implicit in the shell parametrization at lines
   306-314 and is used correctly at lines 714-758; the suggestion only
   makes the new theorem's normalization unambiguous in place.

3. **Lines 513-520: state M>0 in the energy lemma.** Replace the
   hypothesis at line 516 by a phrase such as “若存在 M>0，使其行积分和
   列积分均不超过 M”. The displayed conclusion divides by M; the
   unqualified statement otherwise admits the undefined M=0 case.
   All actual applications use positive M=HM0 or M=4e^2k, so the
   main results and their proofs are unaffected.

These are minor self-containment/formal-domain notes, not hidden
corrections to the accepted mechanism. No other actionable mathematical
change is required by this review.

## 3. Operator, finite cover and actual increments: lines 51-404

The operator is consistently centered, axis-parallel, Lebesgue and uses
one common positive radius in every coordinate. The definition takes
all positive radii and general L1 input. The level-set regularity proof
uses L1 translation continuity and inner regularity without assuming
finite level-set measure in advance.

The increments are defined after each actual selection and frozen; their
nonnegativity makes eligibility irreversible and radii nonincreasing.
The first-rejection argument at lines 192-203 correctly uses only a
prefix whose radii are at least that of the rejected cube, before enlarging
to all selected terms. The factor (beta+1/a)^n/tau has the correct
normalization. Zero increments do not invalidate the argument.

The load normalization chooses the last item among the smallest positive
contributors. All other positive contributors precede it, its captured
sum is one nonnegative subset of its real eligibility test, and its
self-weight is at most a^(-n). The uncaptured radius bound 1<=r_j<t,
common scaling z=y/(a+1), Jacobian and density inequality all have the
correct directions. The coordinate corner lemma handles full-face
equalities at isolated radii and does not call those faces null.

The finite witness integration at lines 362-376 applies to a genuine
finite subcover of an arbitrary compact strict level subset. Thus
selected aQ_j remain true witnesses for the original f. Its scalar
a-to-infinity limit fixes n,t,tau and never assumes histories converge.
The positivity proof includes n=1 and the empty-product convention.
This genuinely establishes Theorem 1 as written, with the full input
and radius scope.

## 4. Original scalar preparation and actual-measure bridge: lines 407-533

The uniform integral estimates and the W_n derivation retain the integral
deduction in the correct direction. The elementary absorption at L>=18
and the stronger L>=64 bounds have the same hypotheses as the locked
source; no numerical optimization was silently promoted to a bound.

At lines 481-510 the measure is built from the same original uncaptured
increments. It is not a separately selected density. Reflection and c>1
give set containment in the common finite corner envelope. The inclusion
D subset H_a is pointwise and yields exact alpha^(-n) scaling on every
Borel subset of D; it is not an unproved measure-convergence step.

The energy identity at lines 524-532 has the right two terms and bound
directions. The second term is controlled by the ambient column degree
using mu<=lambda0; no actual-mu nonconcentration hypothesis is missing.
Disjointness of target and source sets is used exactly when paying only
one total deficit. Apart from stating M>0 explicitly, the lemma is complete.

## 5. Main R042 theorem and its explicit power threshold: lines 535-918

The manuscript actually includes all the new proof, rather than referring
to an internal accepted source:

- **535-665:** the changed x=L-s scale, strong geometric and real-kappa
  guards, fixed-face target statistics, positive target mass, bounded
  Gaussian compensation, coordinate-square error and every deleted
  product are established. The correction near s/4 is bounded, not
  incorrectly treated as o(1). The source corner complement is paid
  by a Taylor power bound; no global source-bad-mass argument is reused.
- **667-697:** the normalized side lengths have a_i>0 and decrease
  with target radius. This yields the pointwise complete-column Y
  envelope, including the factor from deleted products, the full 2n
  face count and interval endpoints. It never differentiates r_-(rho).
- **699-760:** the pinned source coordinate is explicit and its
  squared factor is dropped only upward. The n-1 free-coordinate
  second moment retains the positive cubic term. The paper proves
  the coarse 144 bound directly, including a concrete s>=8 bound for
  its exponent compensation. It therefore does not leave the moment
  merely dependent on the separate private C=48 audit.
- **762-788:** the high-column cutoff removes energy by integrating
  raw columns, then Y, then Y^2/(HM0). The direction is correct and
  does not require a small deleted source mass. The exact H,c_E,C2
  comparison is preserved.
- **790-833:** lambda(A)<=nd is the correct bound for subtracting
  tau*lambda(A). Finite-a energy, degrees and target mass scale by
  alpha^(-2n), alpha^(-n), alpha^(-n), respectively. Each source
  item has strictly larger radius than the same actual target item,
  and the strict center guard puts it into that target's one original
  eligibility test. No independent history or selector is introduced.
- **835-889:** tau is positive and its Taylor-eight power upper bound
  makes the energy lemma applicable. The finite load bound is passed
  through the same original witness/cover interface, so general L1
  and all common positive radii remain included. The new scalar cost
  uses x=L(1-d^2), not the old W_n parameter. The four-term R(s)
  pays log-x change, constant one, A and the square-root cross term.
- **892-918:** base inequalities and decreasing powers cover every
  s>=2^36, hence every log n>=2^108. No finite-box extrapolation is
  used. The abstract and Theorem 2 state the same threshold and
  strict final inequality as the proven scalar remainder.

The new negative correction's order comparison is correct; it retains
leading coefficient 1. It is not a statement about globally optimal
parameters or dominance over every finite-dimensional certificate.

## 6. Older explicit bounds and comparison: lines 919-1400

The manuscript explicitly resets x=L-3log L before the earlier energy
argument. Its row/column products, global source bad-mass bound and
finite-a transfer match the locked R035/R036 chain. The source corner
tail is repaired inline at lines 993-1003 using Ld/32>=2l and
sqrt L>=16l; the withdrawn u^(-2)-only interpretation is not used.
The fourteen envelopes, their true ratio test, zV monotonicity and
exact rational base bounds are all present. The l>=128 theorem is
therefore proved in the text, not just cited as checked externally.

The R034 family at lines 1269-1364 preserves the domain, both positive
factors, shell normalization, exact exponential tail direction,
finite-kappa loss and one actual variable-P argument. Positivity and
qstar slack are present, which justifies its first *strict* inequality.
Its q=1/20 hand-picked example uses a valid rational Taylor certificate.

The min statement is restricted to each formula's valid domain. The
R034/R036 deficit comparison is on the same parameter t=n/(L-3log L).
When comparing the R042 and earlier energy formulas at different t,
lines 1386-1393 retain all A, B, deficit and cross-root terms. The text
explicitly forbids determining the full difference from the deficits
alone. Taking minima with applicable complete family functions and 2^n
is legitimate; no unavailable formula is extended below its threshold.

## 7. Real-history barrier and the fixed-level spike: lines 1402-1520

The paper includes the R037 restoration proof. It conditions on full
shell faces before restricting the cap, retains the M0^2 prefactor,
uses a global uniform slab bound independent of the number of refined
blocks, and assigns constant radius labels above each block's norm.
Each original small cube is exactly a times its grid cell, so real
increments remain disjoint and weights are exactly a^(-n).

The uniform discrete estimate bounds *all* source items, including self
and future sources. This is stronger than each previous-only test and
makes the largest-radius history legal. The positive last minimum-radius
cube is actually restored with common star support and zero captured
sum. The quantifiers fix geometry before a0 and cover every sufficiently
large real a. No relaxation or selector convergence is substituted for
this construction.

The resulting obstruction is explicitly confined to the same algorithm's
universal pointwise W bound in the stated near-linear t window with
tau>=1. It does not claim an operator lower bound or exclude other
joint cover/witness methods. The final spike paragraph fixes n, mass
and lambda before shrinking support and correctly gives both level-set
inclusions. Line 1519 expressly distinguishes this fixed-level ratio
from a whole weak norm; the source's former possible scope confusion
is absent. Jointly witnessing all finite stars remains legitimate
because their open intersection contains a neighborhood of zero.

## 8. Standard results, original results and public scope

Known results A-C have a separate alphabetic theorem counter and direct
citations. The manuscript retains Aldaz's fixed-epsilon quantifier,
uses the generic 2^n bound only as a Lebesgue upper comparison and does
not claim its sharpness for Lebesgue cubes. The paper's own theorems
are separately numbered and proved in the manuscript. Standard tools
such as Jensen, Markov, Chebyshev and Tonelli are used with their valid
measure and positivity assumptions.

There is no claim of global public priority, sharpness, a C_n lower
bound or a smaller leading coefficient. No external theorem is silently
imported from Euclidean-ball estimates. This review uses the previously
audited citation scope and does not undertake a new literature survey.

## 9. First-pass verdict

The mathematical translation is complete and faithful to the accepted
full and threshold proofs. The three minor definition notes in Section 2
are the only requested edits. No production computation, recompilation,
new literature verification or repetition of the original numerical
certificates is required for this verdict.

## 10. Final incremental lock: all notes resolved

The author supplied the original source snapshot at
`/Users/zhengzhihao/Desktop/T/cube-cover-manuscript/cube-dim-7-bounds/upper/round041/manuscript_notes/stage_e523_main.tex`.
Its actual SHA256 was checked to equal the full-reviewed e52394... hash
above. The current 1559-line `paper/main.tex` was independently hashed as
`94c5e61d3bf86b16736e328ffe3e412d46df957dd2b8f8760939edc8831e8692`.
A whole-file Python difflib comparison showed exactly four change hunks,
which constitute only the three requested definition completions:

1. The ell_t and I_n(t) definition is placed before its first lemma;
   the old later `:=` is changed to `=`. The exponent and normalization
   remain exactly those of B_n(t).
2. The energy lemma now explicitly assumes a positive bound M>0.
3. The full source measure is explicitly defined with density
   2^(-n)||z||_infty^(-n) on ||z||_infty/t in 1-I, and
   1-I=[1-3d/4,1-d/4]. Its complete-face radial factor d rho/(2rho)
   and n-1 free uniform coordinates are stated. Closed radial endpoints
   instead of open ones have zero radial measure and change no estimate.

No mathematical inequality, constant, theorem domain, proof branch,
comparison function, bibliography claim or spike-scope sentence changed.
The full mathematical acceptance therefore carries to this final version,
and all three notes are closed. The locked R042 full and threshold
dependencies remain unchanged. No old mathematics, production certificate
or n=7 calculation was rerun, and this auditor did not compile or inspect
the final PDF. The author's compilation/layout report remains separate
from this final mathematical-source acceptance.
