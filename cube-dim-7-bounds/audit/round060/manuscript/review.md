# R060 manuscript: independent full mathematical source review

2026-10-28. **ACCEPT.** The final manuscript contains complete proofs of the R056 explicit operator bound and the R059 fixed-wideband eventual operator bound, with their domains correctly separated:

    log n>=2^20: C_n<=Phi<T_n-n/40;
    all sufficiently large n: C_n<=Phi_w<T_n-n/10.

The second result has no explicit N in this manuscript. No R063 threshold or R061 candidate is imported. The object remains Lebesgue measure, the centred axis-parallel cube maximal operator, arbitrary L1 inputs and all positive common radii. Leading coefficient 1 is unchanged. No blocking mathematical, quantifier, definition or transcription error was found; there are no requested author corrections. The coordinator retains final integration/proof responsibility.

Final locked TeX: `/Users/zhengzhihao/Desktop/T/cube-cover-manuscript/cube-dim-7-bounds/paper/main.tex`, SHA-256 `a5ffdfac3871d4ecd3d7a8e27ecb484019d62f309ce2d7a6e434754e8e17e445`.

## 1. Read scope, sources and execution limits

Every line of the mathematical source was read: abstract/introduction/theorems, all finite-cover/cap preliminaries, the new local d<=1/3 proof, complete explicit R056 chain, complete wideband chain/common eventual conditions, all complete-root comparisons, low-domain R034 proof, real-history R037 barrier/spike qualification, and bibliography. Accepted sources were used for comparison, not as a reason to omit reading their transcribed proofs.

Actual frozen dependencies were checked:

| Source | SHA-256 |
|---|---|
| R056 full proof | `6986ed9ec48b46f3136a9a56722d0092db4f2942d52e35c9f4de37eb279f399a` |
| R056 full audit | `51379eca8cb5a6a7714e89349e96c8c5cfc6c0933ec312d289e039788c10f83f` |
| R059 full proof | `fe95f369634bbc753f5e44a32d6a2b74b288748ebd4bfe1f6c3cac43b086bd8f` |
| R059 full audit | `c8aff39bf68749d25f9978cef11be82ccb7907ec7e72c5678ac66b06a65a33d4` |
| Expanded one-coordinate lemma | `97430af9980c2de5d105bc4ea86667573c5c28e0ff40a92f1c613cd95a243eec` |

Full absolute paths and actual hashes are recorded in `sourcehash.json`. This auditor wrote only the private `audit/round060/manuscript/` directory. No author/main source was edited; no Git, production diagnostic, numerical parameter search, old grid, old local-J baseline, finite-gate rerun, MATLAB or compilation was used. The writer's new 60 checks were not run. PDF compilation/layout is owned by writer/coordinator and is outside this source review.

## 2. Original operator, finite witnesses and actual measure

The strict-level open-set proof handles arbitrary L1 input and possibly infinite level measure correctly. Its compact finite cover uses each point's genuine positive witness radius; selected dilates remain witnesses of the same input. Eligibility sums are nondecreasing, selected radii nonincreasing, and actual increments/weights stay fixed. The unselected-cube argument gives the correct cover factor.

At a positive load point, the terminal contributor is correctly chosen as the last selected item among the minimum positive-contribution radii. Captured predecessors cost tau, the self-weight costs at most a^-n, and other contributing predecessors have normalized radii in [1,t). Their actual disjoint increments give the cap density domination. Coordinate monotonicity of H is proved with the full signed-face formula and strict face counts; exceptional radial equality values have zero radial measure. The finite corner integral, its scalar a-limit and B positivity include n=1 correctly.

The actual measure is the original uncaptured predecessor measure, not a separately sampled shell distribution. On the new supports, lambda_a=alpha^-n lambda holds exactly by D subset H_a. The energy lemma uses a single missing measure and disjoint target/source sets; its expansion pays each missing side once. Its degree denominator M>0 and actual-row hypothesis are explicit. No concentration or full-face probability law is imposed on actual mu.

## 3. Local expanded domain, baseline and full matching

The manuscript proves the local lemma directly for d<=1/3, R>=1, rho>=4r and 0<g<=min(d,R/2). It does not extend the old d<=1/16 theorem by name. The exact Psi/J expression includes a removed tail beyond its onset and, when m>g, the corresponding central deletion. Its lower J bound and Psi partial derivatives justify the numerical derivative constants in the expanded domain.

The capture onset is a genuine kink. The proof uses finite smooth branches, local Lipschitz continuity and absolute continuity, so the derivative is an a.e. bound that can be integrated to every matched q. The pin-crossing join at |w|=R is C1 because Psi_s(R,0)=0. The correct baseline is qb=r-d+h, where R=r; at the matching centre r-d the exact EJ=1 cancellation would not hold. Both constructions use qb.

Baseline free-coordinate Fubini cancellation, the elementary second moment and Jensen/log inequalities give the 3/6/12 constants. The fixed target-face coordinate is deterministic and contributes log in [-g/rho,0], not a free uniform coordinate or a variance term. Individual J precision is independently uniform over coordinate values and matched q; it is not obtained from an n-coordinate failure union or by raising a coarse J enclosure to the nth power.

The two source-sign-count endpoint events retain their n-1 means. Monotonicity of T_N and ell_* gives the lower drift 2h+1/n and the exact reverse-pair difference 2h(1-epsilon)(d-h)/[r(r-2h)]<=2hd. The resulting whole-q normalized count calibrates the pin factor in the full face sum. It does not independently maximize signs at different q. Combined with the product bound and individual J, it produces the stated Ulo/Uhi in their correct directions.

## 4. R056 explicit domain and weighted degrees

For d=3/10 the target interval is trimmed by 2^-20. The finite h guard makes every target's matching interval lie inside the full source band with exact width 2h, including at endpoints. The full source-band G and b are defined before D restriction and b is extended by zero only outside that band. The final cap-restricted kernel, degrees, 0<=K<=1 and disjoint supports are explicit.

The target baseline mean/variance and two count Chebyshev bounds are uniform on full target faces. The corrected original corner bound is per signed face without the invalid averaged half-factor: zero on the negative face and ell_t(r)^(n-1) on a positive one. Four bad events give the stated retained mass 63/64; no independence between them is needed. The four-term log polynomial is explicitly a rational lower mass bound rather than the exact logarithmic mass. Since d=.3, K_w<=n<=4nd is automatic in this construction only.

The source inverse-shell coefficient retains every pinned factor. For a fixed source, a_i=kappa*rho-|z_i|>0 and target capture count divided by 2ell_i is at most one. Strict capture at equality only decreases that count. Thus the full target shell coefficient is <=n product ell_i and yields the new column 2k, rather than an independently deleted-coordinate 4k bound. This is pointwise and survives target-good/cap restriction.

The exact phi first/second weighted moments have the correct double-tail factors. Free coordinates are independent within each fixed captured full source face; the pin has zero variance. The inverse-capture normalizer is J_i. The Q upper polynomial and rate bounds give a common mean gap x/400 and per-face tail 1600/(3x). The proof mixes these **probabilities**, not merely within-face variances. Full-source retention is established before cap deletion; no independence is reused after D or soft conditioning.

The original cap uses the essential exponent x/25, proved by (74/75)*(3/74)=1/25, and pays the gap 1/25-chi=1/200. The tau exponent is separately 49/20000. The Taylor-eight bounds are exactly the cap L^-7 and tau L^(-13/2) envelopes. Target means/variance, count rates, individual J, h/drifts, epsilon/h and inverse-retention upper bounds all have the stated decreasing powers; geometric lower margins increase. The displayed base inequalities therefore cover the full L>=2^20 domain with no unquantified eventual condition remaining.

The upper row is below (21/8)k, columns <=2k, cap loss and actual threshold each cost e*k, and the reserved epsD gap is distinct. The exact gamma matches the main theorem. The scalar price 7/93, A term and positive square-root cross term are all retained; the Taylor-eight P(L) powers L^-8,L^(-15/2),L^-8 decrease. Its strict base slack is the accepted exact rational value. Thus the main operator domain is a genuine complete threshold, not a scalar-sorting start.

## 5. Wideband geometry, five events and active-source restriction

The manuscript changes the complete fixed tuple and explicitly reuses local notation. It does not transfer the preceding L0. On the full wide source band R>=R0-h>1 eventually, so b/G are defined even for inactive sources. A source with any matched target has the stronger R>=R0+delta-2h>=R0; that bound is used only for active-source pin moments. Inactive columns are zero. The inset ensures full matching width 2h and source/target radial separation. The local lemma applies on every active matching interval; baseline remains qb.

Here d<1/4, so K_w is **not** automatic. The paper includes it as the fifth bad event, proves its mean/variance and the 1/(4nd^2) bound, and allocates five uniform o(1) failures by epsA/5. Baseline, both count endpoints and target cap have the correct full-face laws. Their union gives the exact full target annulus mass and the strict lower m0. This is an ambient statement only, with no law assumed for actual increments.

The positive-log S64 and geometric tail formulas are correct: the first omitted exponent is 129. The uA and uQ substitutions represent log(rhi/rlo) and -log(1-d/(2R0)), respectively. The asserted finite rational inequalities give the mass/pin bounds and moment rates. Independent new six-term positive-log certificates, with their own rational tails, also verify these two strict margins; this does not rerun the writer's 64-term certificate.

On active sources the pin moment uses R>=R0 and R/q<=R0/qmin. The per-face upper mean 1/9+(293/20000)x and variance x/900 give gap x/40000 once x>=40000/9. Per-face Chebyshev followed by mixing is again valid; no total-variance shortcut is introduced. The complete product, individual J, normalized count and column-two formulas transfer under their explicitly reverified domains.

The wide source mass is paid with Csrc=(qmax-qmin)/qmin. This coefficient is larger than d for this tuple; silently importing the old nd coarse mass bound would be unjustified. The manuscript correctly uses Csrc and retains the n-1 factor through [(1-1/n)eta-chi]. Its exact eta-chi>1/100000, followed by eta/n<=half that gap, supplies a positive source-cap exponential. The separate tau exponent is 741/1000000>0. Each exponential payment is a fixed polynomial times exp(-cL), so both can enter a common eventual domain.

The complete common-guard list includes geometry, baseline, K_w, both endpoint counts, target cap, individual J, calibration drifts, inverse retention, source cap and tau, all uniform in radii and faces. Their limits give one absolute integer start for this fixed tuple; the final scalar error can be included in the same finite list. No explicit start is computed or borrowed from R056/R063. The chosen gamma0 lies strictly below (1-e)Gamma0, and the complete price 3/97 leaves the exact positive scalar slack 39/1212500. Thus the eventual n/10 claim is justified, without claiming the c_raw endpoint.

## 6. Same actual history and all-scale operator passage

For both tuples actual target items are below radius 3 and actual source items above 4, so all included sources precede the same actual target. The displacement guard for every original dilation mathfrak a>=8n^4 is n^-3-(beta+2)/(mathfrak a+1)>0. An included source increment piece costs at most its original weight; bounded b<=1 therefore preserves that target's eligibility tau, without an inverse-theta factor.

The exact lambda_a scaling multiplies energy by c_a^2 and degree/target mass by c_a, c_a>=1. Applying the one-missing-measure energy identity cannot reduce the resulting absolute deficit. Both full root budgets therefore apply uniformly to all original finite histories at each permitted n. Holding n,t,tau fixed, only scalar inequalities are passed to the dilation limit; no selector/history convergence, independent history splicing or altered input is used.

The original compact strict-level finite witness cover and Lebesgue inner regularity apply to every general L1 input and all positive common radii. There is no hidden bounded support, radius lattice or scale cutoff. The wideband section repeats this interface rather than inferring an operator claim merely from a model objective or accepted numerical certificate.

## 7. Complete cost comparisons, R034 and real-history barrier

The two older exact formula costs are explicitly defined and positive on their table domains. Their lower envelopes retain the complete B-minus-deficit component, so the new explicit Phi<T-n/40 is strictly below them on those domains. This does not require a removed old theorem to remain an unproved self-claim: C_n<=Phi<Phi_j is proved directly.

For Phi_w versus Phi, the B lower envelope correctly retains -log(1-7/100)>=7/100. Thus Phi>T-n(gamma-7/100+4/L); its bracket is below 1/10 already at L=2^20 and decreases. The eventual wide upper Phi_w<T-n/10 then orders the **complete** root costs on their common domain. Both A/cross terms and the positive parameter price have the correct directions. This does not assign an explicit start to wideband.

The low-domain R034 theorem/proof is intact and self-contained: negative-face target mass, shell overlap and face count, corner subtraction, kappa-shrink product and availability losses, same-target variable capture set, exact finite-a scaling and an interior q_* to obtain strictness are all present. Its Y=64 q=1/20 instance pays the exponential by a positive finite Taylor sum. W_n has its own old t and is not used as the new-t scalar cost.

The retained R037 barrier constructs real histories with correct labelled radii, compact gaps, uniform boundary layers, full incoming grid load and a positive final minimum-radius item. It constrains pointwise load bounds for this selection algorithm in a fixed nearlinear t-window. The text explicitly excludes a C_n lower bound. Its common-spike ratio tends to one at fixed n,mass,level; it expressly does not claim the entire weak norm tends to one. The distinction from joint input/coverage information remains intact.

Known results A/B/C remain distinct from the original numbered theorems. No public priority, global parameter optimum, R058 ideal value or complete numerical ranking is a mathematical premise of this manuscript.

## 8. Final formatting delta, actual check and outstanding work

The initially read full source had SHA `93b9af22489aa2d8fc55d0ab3d01373c89175aa60efd88996a35d60268215aa4`. The writer then changed only the separator after tqmin>=4rhi in the wide common-guard display from a quad to a mathematical line break. The exact patch has SHA `035dbd27fac48df2c505f35e497245c2cd52d48403922105c5eed86bffb745c5`. This auditor inspected it and independently reversed that single byte-level replacement in memory, recovering the original 93b9 hash exactly. Thus the final a5ff identity changes no mathematical content.

`source_check.py` was actually run locally, PID 92655. Seven new compressed Fraction/source certificates PASS: six-term positive-log mass/pin enclosures with rational tails, the fixed wide budget, the source-width distinction, exact net slack, and both full-root comparison margins. The same run verifies 93 labels and 58 references, with no undefined/duplicate labels or undefined citations, and records all actual source hashes. A preliminary read-only probe also checked the six-term enclosures and reference inventory before this registered combined run. No author 60-check or old mathematical baseline was executed.

The final full ACCEPT is locked to `a5ffdfac3871d4ecd3d7a8e27ecb484019d62f309ce2d7a6e434754e8e17e445`. There is no outstanding mathematical or definition correction. Writer/coordinator remain responsible for actual compilation/rendering, and later explicit-wideband developments remain outside this frozen manuscript review.
