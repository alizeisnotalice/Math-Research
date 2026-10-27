# R024 independent pure audit: unique minimum of the proved family

2026-10-27. **Accepted; no source mathematical erratum found.** The all-a>1single-minimum theorem is correct for the specified two-parameter upper-bound family. The fixed-a threshold, derivative ratio, sign partitions and finite root-localization contract are valid. This audit has not evaluated F36, endpoint signs, any midpoint, or a minimum value. Parent retains numerical-certificate and overall proof responsibility.

## 1. Frozen sources and inherited theorem

The current proof is SHA25672d7af4ed81f8b456dc8d636ac0251c63f175e05405b5f8581b9f890cd1d3cc3; parameters SHA2561257cfc42eb7e05ea29e59720f0b79ea0f454b1e916b31c68703a36645aae5cd; registration SHA2563d1b66dd0754dcf1a731239600944e5b36ae0ddc97695598c0b824ad9c212707. Their complete text, sources and receipt were read. The accepted R023parameter theorem was also read in full; its allL1/allpositive-common-radii centered-cube conclusion is the inherited proof dependency. This round changes the deterministic parameters, not that operator or its input class.

The inherited formula is C7<=K(a,tau) for every a>1,tau>0, with

    K=(1+A/tau)(tau+B),
    A=(1+2/a)^7, B=(1+1/a)^7(1+7log a).

The new independent clone is atb135609efca22428ac8f6a9c9485d114a4a03846 and writes only audit/round024/parameter_theory. No source/older artifact was changed. Original R023numerical result/manifest bytes are old comparison evidence and are hash-checked without rerunning a row. The initial path query before source creation was unavailable; after the author supplied the frozen source, its full text and hashes were used. No mathematical execution occurred during either stage.

## 2. Exact tau optimization and the one fixed expression

For each fixed a>1, A,B>0 and

    K=A+B+tau+AB/tau.

The last two terms are uniquely minimized at tau=sqrt(AB), giving

    F(a)=(sqrtA+sqrtB)^2.

This is an actual minimum over all positive tau, not a stationary-point heuristic. The chosen threshold is fixed as a function of a and dimension only, independent of f,lambda,radii or increments. The parent theorem is uniform in these variables, so choosing that real threshold is legitimate.

At a36 the source expression is exactly

    F36=[(19/18)^(7/2)+(37/36)^(7/2)sqrt(1+7log36)]^2,
    tau36=(703/648)^(7/2)sqrt(1+7log36).

The rational product703/648is(19/18)(37/36), so the threshold scale is correct. This expression is no larger than the accepted same-a,tau7upper expression by exact minimization. A strict numerical separation against a saved old interval remains a separate test; this audit does not claim its outcome.

## 3. Independent derivative factorization

Let H=sqrtA+sqrtB, L=1+7log a and N=a-7log a. Direct differentiation gives

    (sqrtA)'=-(7/a^2)(1+2/a)^(5/2),
    (sqrtB)'=(7/(2a^2))(1+1/a)^(5/2) N/sqrtL.

Factoring the positive first magnitude yields

    H'=(7/a^2)(1+2/a)^(5/2)(R-1),
    R=1/2*((a+1)/(a+2))^(5/2)*N/sqrtL.

This is algebraically the source P(a)factor and confirms the essential5/2exponent. Because H>0, F'=2HH' has exactly the same sign. The positive square-root/log domains are automatic on a>1. No numerical derivative or symbolic software was used.

## 4. All-domain sign proof and unique global minimum

For1<a<=7, N'=1-7/a<=0 and N(1)=1, so N<=1. If N<=0 then R<=0. If N>0, the ratio factor is<1 and L>1, so R<1/2. Thus H'<0 throughout this entire first range. No claim that N itself stays positive is needed.

For a>7, N is strictly increasing and N(7)=7(1-log7)<0. One elementary independent proof log7>1 uses integral_1^2 dt/t>1/2and integral_2^4 dt/t>1/2, obtaining log4>1, then7>4. The source's alternative e<3<7argument is also valid. Since log a/a tends to0, N tends to infinity and has a unique zero a0>7. Wherever7<a<=a0, R<=0 and H' remains negative.

On a>a0, both N and L are positive. The derivative of S=N/sqrtL has the sign of

    2(a-7)L-7N.

Here log a>1 implies N<a-7, while L>8. Therefore

    2(a-7)L-7N >(2L-7)(a-7)>0.

So S is strictly increasing wherever it is positive. The other positive factor((a+1)/(a+2))^(5/2)also strictly increases; their product R is strictly increasing on this range. The multiplication argument is deliberately restricted to positive S. It would be invalid to infer the product's monotonicity just from increasing factors on the earlier negative-S range.

The source additionally proves the stronger S'>0for every a>=7 by its displayed T>=7(log a-1)/2; expanding T and bounding7/(2a)<=1/2and7/a<=1 verifies that estimate. It is consistent with the independent restricted argument above.

Finally R tends to infinity: the elementary bound log a<=2(sqrt a-1) implies N>=a/2eventually and sqrtL is at most a constant times a^(1/4), while the ratio factor tends to1. At a0, R=0. Continuity and strict increase after a0therefore give exactly one a_*>a0with R=1. The derivative is negative on(1,a_*)and positive on(a_*,infinity). Consequently F strictly decreases then strictly increases and a_* is its unique global minimum over **all a>1**. Combining with the unique tau minimum gives exactly one jointly minimizing pair within this family.

This proves existence and uniqueness without a numerical root bracket, optimizer or sample grid. An unexamined boundary a→1cannot hide a better infimum because F is strictly decreasing throughout the interval from1to a_*.

## 5. Fixed endpoint/bisection contract

The proposed signs R35-1<0<R37-1are prerequisites, not analytically or numerically certified by this audit. If the separately authorized outward arithmetic certifies both, the unique-root theorem places a_* in(35,37). Each midpoint is an exact rational and receives one strictly separated sign enclosure; replacing only the endpoint with matching sign preserves the root bracket. A zero-containing interval is unresolved and must stop the stage. It cannot be assigned by a midpoint heuristic, interpreted as exact equality or trigger an automatic precision/bracket change.

After80successful halvings, width is2/2^80=2^-79. The final exact midpoint a_c is within2^-80of a_*. Evaluating F(a_c)gives a legal direct C7upper expression simply because a_c>1; its legality does not depend on a small root error. The interval for this **candidate value** need not contain F(a_*), so it is not a minimum-value enclosure. A separate interval/Taylor argument would be needed for that different claim and is explicitly outside the contract.

The prescribed84function evaluations are exactly F36, two initial Rsigns, eighty midpoint Rsigns and F(a_c). Cached intermediates do not create new parameter rows; actual library work must still be reported honestly. The positive-root implementations r^3sqrt(r)and r^2sqrt(r)correctly implement7/2and5/2powers at positive rational r. Every rational input is formed exactly before conversion to the ball library.

The strict rational grid upper rule floor(U*10^12)+1over10^12is valid even when the outward upper endpoint lands exactly on the grid; it remains strictly larger than U. An old-grid strict improvement may be reported only if a new upper endpoint lies below the saved old lower endpoint. A nonseparation is not proof of equality or failed validity and does not authorize more rows.

The main has separately released the numerical owner; this audit neither starts that task nor treats its future endpoint/step outcomes as already passed. Partial failure should preserve honest data and stop as specified. No endpoint or midpoint expression has been evaluated here.

## 6. Scope of the optimum and final acceptance

The inherited estimate holds for every deterministic parameter pair uniformly over all generalL1inputs and all positive common centered-cube radii. Choosing an exact a_c,tau_*(a_c), or taking the infimum of these proved upper expressions, preserves the same weak-type operator and hypotheses. No Green-floor assumption, spatial restriction or additional factor is introduced.

The unique minimum is **only the minimum of this proved expression family**. It is not the unknown best C7, not an optimum over all covering arguments, and not a C7lower bound. A lower enclosure for F(a_*)would still be merely a lower bound on these upper-certificate expressions. The old R023grid and any new legal candidate are separate upper certificates, not operator counterexamples or literature-optimality claims.

Strictly accepted: exact threshold minimization; a36closed form; derivative factorization; all-a>1unique minimum of this family; and the conditional finite root-localization/candidate-value protocol. Not certified here: the initial bracket signs, eighty successful steps, a numeric optimum location/value, strict improvement over the old saved interval, or any new numerical C7constant. No source erratum was found. Source hashes, independent receipt and manifest accompany the review.

Only proof/source reading, independent mathematical reasoning and metadata hashing were performed. No mathematical numerical or symbolic software, scalar, bisection, optimizer, sample, phase/DD, MATLAB, commit or push was used. This audit is complete and stops; no additional function or task is introduced.
