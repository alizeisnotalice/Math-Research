# N3 round 003 — three layers and staircase capture

2026-10-25. **No new three-dimensional bound.** The optimized three-layer finite union is strictly below the inherited accurate R3. There is a separate rigorous structural finding: keeping layerwise staircase captures increases one stated finite witness union by more than 0.12185552575179436. That improvement is within a weaker construction; it is not an improvement of C3's best lower bound.

## Mechanism and saved exact constructions

Background is Lebesgue; all witnesses use centered axis-parallel cubes with one common radius. Start with a=(1+u/2,(2+u^2 v)/4,(1+u)/2), u^3=2, v^3=3, directly from the original certificate, and periods Pbase=(3a1,a2,a3). Three orbits have positions P*(k+phi_b), k∈Z^3, b=0,1,2. During discovery, phi_0=(0,0,0), phi_1=(1/3,z3,z4), phi_2=(2/3,z5,z6); P=Pbase*exp(z0,z1,z2), masses (1+z7,1+z8,1-z7-z8), total mass 3 per orbit-index triple. The independent transverse phases, mass redistribution and periods are optimized together. This is a new three-cell family, not the old two-orbit 143-key optimization.

A capture stores a separate inclusive integer rectangular block [l_b,h_b] for each orbit. Different transverse block lengths and starts across b give genuine staircase shapes in original layer coordinate (3k1+b,k2,k3). Every selected atom in every orbit enters the same cube. This follows the actual missing-subset mechanism in the four-dimensional report rather than assuming the whole selected set is one original-index rectangle.

`round003_certificate_input.json` saves exact rational parameters and every selected block key. Each key has 18 integers: nine lower endpoints orbit by orbit, then nine upper endpoints. Empty orbits use l=(0,0,0), h=(-1,-1,-1). The optimized survivor has exact periods (4.921704158643,1.072321649545,1.129926430968); its saved phases and masses are terminating decimal/rational values with total mass exactly 3. Discovery's 1/3 and 2/3 longitudinal phases were rounded to the saved 12-decimal rationals in this certification; that tiny geometry change is included in every bound. No parameter was derived from decimal digits of R3.

For each key, M=sum_b w_b product_j max(0,h_bj-l_bj+1), r=M^(1/3)/2. The receiver box in the full cell ∏[0,P_j] is the coordinate intersection of [max_b P_j(h_bj+phi_bj)-r, min_b P_j(l_bj+phi_bj)+r]. Only nonempty orbits enter extrema. All receiver boxes are retained in their full three-dimensional joint union; neither marginal multiplication nor a reflection reduction is used.

## Executed finite batch and recapture

The initial 37 structured cases comprise baseline plus monotone, opposite, and rotated transverse layer patterns, amplitudes .01,.04,.1,.2, and mass perturbations (0,0),(.1,-.05),(-.1,.1). Each uses 4096 scrambled Sobol receivers, seed 20261025303. At each receiver all atom-entry radii <=2.8 are tested, with cube counts computed separately in each orbit. The ±3 entry-index window covers this radius range for all registered periods. There is no claim about radii >2.8. The 37 cases generated a pooled 225 integer keys. Baseline joint union was 1.7939332277947664, closing against the original algebraic formula numerically.

The four highest nonbaseline starts each took 350 full-union Nelder–Mead evaluations. Each survivor was independently recaptured at 16384 Sobol receivers, seed 20261025304. New-key counts were 22,24,22,31. Their joint unions were respectively about 1.79361138227957, 1.79344021223088, 1.79338656803006, 1.78702181102913. These extra keys did not visibly increase the respective unions beyond floating-point last-digit noise; this is not an exclusion of omitted sets.

The best survivor then took the registered last 400 full-union evaluations, followed by a final 16384-point recapture, seed 20261025305. That last recapture added three keys. The final active family has 232 witnesses and floating value 1.7938623489136238 before rationalization. All five optimization runs stopped at their evaluation cap, so no local optimality is asserted.

A serialization error occurred after the original 1800 optimizations and final capture had finished: NumPy int64 keys could not be written as JSON. Four completed 350-evaluation records survived. The pre-registered recovery kept those four records, converted keys to builtin integers and replayed only the same deterministic final 400-evaluation path. Actual optimization execution was 2200 evaluations; the unique optimization path remains 1800. The registered diagnostic replay reconstructed the same 37 case/key families and compared rectangular versus staircase subfamilies. Original and replay captures therefore total 74 initial-case captures and seven follow-up captures, including the two repeated follow-ups needed for recovery. No new optimization start, parameter mechanism or further batch was added.

The first script's optimization guard had also assigned zero to three monotone .2 scan cases whose second phase is .4. The diagnostic replay corrected that guard to .5 and computed those registered cases' actual unions; it did not alter the four selected .01 optimization starts. The best structural diagnostic below comes from one of those now correctly evaluated larger offsets. These implementation corrections and recovery were registered before their execution. The diagnostic/recovery pass took about 12.8 s; the original run's wall time was not preserved by its failed output serialization.

## Accurate results and a concrete staircase

`round003_certify.py` reads only the saved exact input. Fraction arithmetic brackets every cube root at denominator 10^24 and verifies its cubed integer inequalities. Inner and outer boxes give nested full unions, computed by exact recursive coordinate decomposition, divided by total mass 3. Each upper endpoint is directly compared with the **original R3 rational lower endpoint**. The decimal enclosures below are outward displays; full exact fractions are in `round003_certificate_result.json`.

|Saved finite union|Keys|Certified ratio enclosure, outward display|Comparison|
|---|---:|---|---|
|Optimized three-layer survivor|232|[1.7938623489139913956485414,1.7938623489139913956485444]|upper < original R3 lower|
|Selected staircase diagnostic|99|[1.7035932629883014281370510,1.7035932629883014281370555]|upper < original R3 lower|
|Its selected original-index rectangles|53|[1.5817377372365070627814893,1.5817377372365070627814931]|upper < original R3 lower|

The optimized upper endpoint falls below the original lower endpoint by more than 0.00007087888077928. The staircase diagnostic uses transverse phases (.2,0),(.4,0) and masses (1.1,.95,.95). Splitting its captured keys by whether their selected count fills their minimum enclosing original-index rectangle gives 53 rectangular and 46 nonrectangular keys. Its exact expanded-union lower endpoint minus rectangular-union upper endpoint is strictly positive and greater than **0.12185552575179436**. This compares two stated finite unions at the same rational input; the 53 keys are not claimed to enumerate every possible rectangular witness.

One explicit six-atom staircase in `round003_staircase_diagnostic.json` is:

- orbit 0: k1=0, k2=0..1, k3=0..1 (four atoms);
- orbit 1: k1=0, k2=0, k3=0..1 (two atoms);
- orbit 2: empty.

Its exact mass is 4(1.1)+2(.95)=6.3. In original layer coordinates its minimum rectangle has eight points, so the subset is genuinely nonrectangular. Its spatial spans are approximately (a1,a2,a3); adding the missing second-orbit row makes the transverse span 1.2a2. This explicitly explains the geometric loss hidden by a one-rectangle selection rule. The whole 99-key certificate, including this example, uses the saved rationalized positions.

## Strict levels, finite mass and explicit index margin

For all three saved unions, active orbit indices have lower/upper bounds within [-2,2]^3. The optimized and staircase families have index_min=(-1,-1,-1), index_max=(1,2,2); rectangles have index_min=(0,-1,-1), index_max=(1,2,2). The verifier computes **b=2** in each case from all keys and asserts every active block endpoint is inside [-b,b]^3.

Receive physical cells P*(q+[0,1]^3), q∈[-N,N]^3. Translation of a witness adds q to every orbit label k. Retain all three original orbit representatives at every label k∈[-N-2,N+2]^3. Their mass per label is exactly 3, including representatives whose negative transverse phase places them outside the named physical cell. Thus finite input mass is exactly 3(2N+5)^3 and receiver union volume is (2N+1)^3 U_tau. No reindexing by physical cell membership is used. The actual verifier checks endpoints at N=0,1,7 for all active orbit blocks (1467,510,207 checks respectively).

To obtain strict superlevels and actual L1 inputs, contract each retained receiver box by tau on every face and reduce its critical radius from r to r-tau/2. The captured atoms were at distances at most r-tau; smoothing them by same-mass nonnegative bumps with infinity-radius <tau/2 preserves containment, while (2r-tau)^3<M. Use one sufficiently small tau for finitely many boxes. The finite weak ratio is at least (U_tau/3)*[(2N+1)/(2N+5)]^3. First N→infinity, then tau→0 gives C3>=U/3. All saved inner-root bounds happen to be strict here, but the contraction argument handles critical equalities as well. This proves valid finite-witness lower bounds for the all-scale, general-L1 definition; none exceeds the inherited R3.

## Actual local checks, proof status, next mechanism

Python executed the registered scan, optimization, recaptures, recovery and exact certification. An independent integer-grid union check with seed 20261025306 compared the recursive union volume to direct unit-cell occupancy in 32 cases, 1..7 three-dimensional integer boxes each; all passed. This tests finite union arithmetic on those cases, not every production box independently.

MATLAB actually ran `round003_matlab.m`, independently formed all boxes and evaluated their full joint unions using a three-dimensional difference grid. It returned optimized=1.7938623489139902 (232 boxes), staircase=1.7035932629883019 (99), rectangular=1.5817377372365067 (53), mass 3 in each. All variables use n3r3 prefixes; no clear, path or shared-working-directory change. These are floating cross-checks; exact inequalities come from Fraction.

- Numerical candidates: five capped optimization survivors, finite receiver sampling and discovered capture keys. No convergence, exhaustive capture or family extremum is asserted.
- Local sufficient conditions: selected orbit blocks, explicit common-radius containment intervals, uniform contraction and finite-index transfer.
- Complete certification: exact lower/upper enclosures and original-R3 comparisons for all three saved finite unions; strict positive staircase enrichment at one fixed input; b=2 and exact finite masses.
- Excluded: only the claim that these stated finite unions improve R3. Their complete periodic superlevel sets and the three-layer family remain open. No new upper bound or whole-family exclusion is claimed.

A distinct next mechanism is an original unit-mass lattice plus a genuinely small free satellite orbit epsilon*delta_(z+aZ^3), optimizing a,z,epsilon. One could first locate z from first-order exposed-face displacement: growth of receiver volume per epsilon must exceed R3 to compensate for mass 1+epsilon, whereas independent new singleton boxes in old holes contribute derivative 1. This differs from nearly equal layer-mass redistribution and was **not run in this batch**.

## Reproduction and integration

Use the supplied external venv Python. `round003_search.py` reproduces the substantive 37-case plus 1800-evaluation search with the serialization fix. `round003_completed_runs.json` preserves the four actually completed first-stage records. `round003_diagnostic.py` reproduces the diagnostic/recovery from them and regenerates the final exact input; `round003_certify.py` independently of search reads the saved exact input and regenerates exact certificates/checks. Run `round003_matlab.m` for the independent floating union check. The scripts import only the already archived local union/root helpers from `witness.py` and `certify.py`; hashes pin those dependencies.

Integration whitelist is precisely the new `round003_` files listed in `round003_source_manifest.json`: registration, report, three Python scripts, MATLAB script/output, completed runs, compact results, staircase diagnostic, exact input/result and source manifest. The optional `round003_keys.json` convenience duplicate is omitted; the exact input contains those keys. Temporary logs and Python caches are omitted. No previous evidence, STATE, total ledger or other project was modified. No commit or push was made. No subsequent experiment is running.
