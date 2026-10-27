# C149: finite NEW canonical density compiler pilot, prepared only

This packet scales the exact C147 pushforward construction to at most100 NEW canonical source patterns. It does not calculate k_h, a pair integral, beta2, M or C6. Existing C147 science and the C146 moment certificate are read-only and are not replayed. Sol supplies the mathematical scheme/code; Terra operates all preparation and scientific runs after main authorization.

## Mathematical law and canonical identities

The source law is exactly C130: E_j iid Exp(1), (A,B)=sum_(j=1)^6 (b_j,c_j)E_j. Its joint Laplace transform is product_j(1+b_j s+c_j t)^(-1). Its rank2 positive density and finite piecewise exponential-polynomial construction are proved self-contained in the frozen absolute C147 PROOF reference in the source manifest. Initializing with two independent columns gives det^-1 exp(-lambda.z) on their cone. Adding each remaining nonzero column v convolves with its positive exponential ray law.

In each existing sector the ray parameter has finite linear endpoints. Comparisons between two nonzero endpoint ratios annihilate v; comparisons with0 or a zero denominator use old boundary rays. Consequently adjoining only the new column direction suffices; no arbitrary fan explosion is assumed away. With m active columns and d distinct directions the final fan has at most d-1<=m-1 sectors, and polynomial degree at most m-2<=4. Eta=lambda.v-1 equals0 is integrated as a polynomial, without singular partial fractions. Repeated columns/directions and zero columns are therefore legal. Zero columns integrate out with mass1. All selected source laws are rank2; singular rank1/rank0 remain covered by the old proof but are outside this pilot.

Canonicalization preserves exact magnitudes. It sorts the six exact rational column vectors and chooses the lesser of that tuple and its simultaneous A/B-coordinate swap. Column permutation is valid because E_j are iid. Output swap is valid for the symmetric k_h(A)k_h(B) integral; for a nonsymmetric future integrand the stored swap map must be applied. No direction-only scaling, unequal-rate identification, independence assumption or change of source weights is allowed.

Terra's frozen9961-pattern table represents all10416 original cone/body pairs. Preparation verifies canonical identities, source membership and exact positive alpha sums. It locates the three C147 passed canonical keys from their original source rows and exact saved-result references and excludes them from NEW compilation. Their laws are referenced through permutation/swap only; neither their density nor their old k values are regenerated. Selected NEW patterns retain all original source lists and their full rational aggregate cone weight.

## Finite selection and coverage contract

Prepare-only selects exactly100 distinct NEW keys. Deterministic strata are (active m, distinct directions d, presence of repeated nonzero columns), with at least one from each available stratum before round-robin filling. The packet must cover m2..6, maximum6 directions, zero columns, repeated directions and identical nonzero columns. A NEW m3/d2 case with one duplicate nonzero column is reserved as an eta0 witness. Its initial independent pair contains one repeated column and the other direction; convolving the remaining identical column has lambda.v=1. The actual scientific run requires a positive resonance-branch count for this witness rather than trusting a metadata label.

No configuration has been produced or density compiled by Sol. The code-only bundle must first be copied unchanged into Terra's designated output directory. One authorized Terra preparation run produces the exact case packet, frozen source manifest and registration. Main and an independent reviewer then inspect those concrete hashes before releasing science. Preparation is Fraction/hash/source bookkeeping only: it imports no SymPy, generates no density and evaluates no old functional.

## Exact density and Laplace certificates

On each sector maintain the complete signed sum P(z)exp(-lambda.z). Its positivity is inherited from the convolution law, not from each summand. Eta!=0 uses the finite integration-by-parts polynomial; eta0 uses monomial integration. Algebraically identical exponents are merged and zero polynomials removed exactly.

For every generated stage, use z=u alpha+w beta on each sector, expand the polynomial in alpha,beta and integrate each monomial for sufficiently large positive Laplace variables. The finite rational sum is canceled as a whole and checked exactly against the product of source factors for that stage. Difference numerator0 and mass1 after whole-expression simplification are mandatory. This handles signed terms and removable poles correctly, including resonance. Numerical density samples and clipping signed terms are not substitutes.

Each passed stage gets its own immutable JSON checkpoint containing all density terms, source-active columns, fan, rational transform, exact zero-difference assertion, normalization and counters. Each case has a separate intent, start, result and receipt. A later independent reviewer can reconstruct the exact rational certificate from those files. Future interval integration must enclose the full density sum before intersecting with its proved nonnegative range.

## Honest finite counters and bounds

C147's mixed counter is not reused. Here a literal endpoint counter increments immediately before each actual upper/lower substitution; both resonant and nonresonant branches perform two. Successful cases assert

    literal_endpoint_substitutions = 2*convolution_term_branches,
    convolution_term_branches = resonant_branches+nonresonant_branches.

Stage attempts and completed exact Laplace checks are separate. A failure retains the mutable counters as observed at failure, so no partial work is silently discarded. A killed/broken worker has explicitly unknown counts, not invented zeros; conservative per-submitted-case budgets still bound the finite source work.

Starting at T2=1 term, a conservative bound is T_k<=2(k-1)T_(k-1): T3<=4,T4<=24,T5<=192,T6<=1920. The same recurrence bounds the endpoint substitutions generated at each addition, so at most4+24+192+1920=2140 literal substitutions per full case suffice before exact term merging. The operational cap2200 allows this bound. No stage may retain more than2000 nonzero terms. A bivariate polynomial of degree4 has at most15 monomials, giving at most28800 polynomial/Laplace monomials per stage; the cap is30000.

Additional operational caps per case are5 stage checks,1,000,000 nodes per explicitly inspected expression,20,000,000 inspected node visits, and8192 bits per stored rational coefficient. These last caps are checked at visible polynomial/Laplace boundaries; they are not a claim that SymPy's private intermediate simplifier graph has an instrumented hard node limit. Exceeding any operational cap is a preserved finite failure and excludes no valid density law or mathematical family.

With at most100 submitted NEW cases the registered ceilings are500 stage checks and220000 literal substitutions. Main chooses1..15 worker processes. Each worker imports its own exact algebra runtime; no shared symbolic context is mutated. One-thread environment and nice>=10 apply. There is no time/CPU ceiling, no compiler, no native kernel and no automatic retry.

## Stop, partial result and cost scope

At most the released number of cases is in flight. A shared stop flag prevents new guarded case launches after any case failure; already running cases finish their finite work and are reaped. Each successful stage/case remains independently usable. Per-case checkpoints are never overwritten. A process crash makes its missing work counts explicit; other checkpoint evidence is preserved. No failure launches a precision, term, pattern or worker expansion.

Pattern weights mean only the coefficients for a future C130 pair integral. A passed density contributes no computed number to beta2. A partial batch reports both selected input weight and passed density weight; it does not claim the selected total was certified. The three reused keys are reference-only. No prefix probability is renormalized, no omitted integration tail is paid by this symbolic compiler, and no actual M upper or C6 statement is produced. The full9961-key compilation and any k/infinite-tail integration are outside the release.
