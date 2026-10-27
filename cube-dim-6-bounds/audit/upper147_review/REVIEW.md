# C147 independent proof/source preflight

Verdict: PASS for main to release exactly the two frozen m3/m6 symbolic density cases. This review does not release any kernel integration, whole-source batch, native query, numerical bound or old six-dimensional replay. The reviewed generator was not launched.

Reviewed source directory: `/Users/zhengzhihao/Desktop/T/Math-Research-cube6-sol-localmax-audit/cube-dim-6-bounds/upper/round147/pair_pushforward`.

## Source pins

The sealed source manifest is 2044 bytes with SHA256 `2500cf19538335300f9c4f2db74c914620acba279af8e182a1f43c96eb9de69f`. Every listed file hash and byte count, the four original C130 source hashes, and all cross-pins in the next packet were independently checked.

- `generate_next.py`: `696e6bdaf629fdfc6404e76c525149fc101b798a83878b9dae43fcc6ba4c0daf`.
- `PROOF.md`: `2ecb770367e181fa153c7cfc1d2bcb7c67f7babc17b6ed7b24b7f4c78adeb22f`.
- `next_cases.json`: `4cb15948d91f799c7cc567b771f0a647047d9f52627dc9d66d6b823400ada7b3`.
- `next_registration.json`: `1b97c07081651eb5f4b6e275403dab185a4321c59f8ec611ef760dfe8a96b0c3`.
- `next_release_requirements.json`: `b17d6dd924b9c0250c7ba170c6fdcb03e41bf48a186f0a8b2fa235662a3916f7`.

The fixed rows exactly match C130 cone0/body(9,16), with three nonzero columns in three directions, and cone0/body(0,1), with six nonzero columns in two directions. Both have rank2, nonnegative rational coefficients and original cone weight `1/171992678400`. Column lengths remain unchanged.

## Mathematical audit

The product Laplace law follows from independent input exponentials. The outputs need not be independent. A rank2 pivot has the stated inverse-Jacobian density on its cone. Convolution with each remaining exponential ray preserves nonnegativity, unit mass by Tonelli, and absolute continuity. These properties concern the complete density, never each signed summand.

For an old sector, the two determinant inequalities supply the correct lower/upper bounds on the ray parameter; zero denominators are tested without division. Every nonzero nonnegative ray has at least one finite upper bound, since the sector lies in the pointed nonnegative quadrant. Empty and zero-length intervals contribute zero. Differences of nonzero endpoint forms annihilate the added ray; their only nontrivial zero line is its span. Original boundary comparisons handle the candidate zero and zero-denominator tests. Hence the old directions plus the added ray suffice, including axes and parallel/repeated directions. Interior samples determine comparisons throughout each open refined sector; boundary overlaps have zero planar measure.

The eta nonzero antiderivative is correct. Its endpoint substitution does not increase total polynomial degree. Eta zero is integrated as a polynomial and increases degree by at most one. Starting at degree zero gives degree at most m-2, including repeated columns/rates. Removing zero columns integrates their probability mass to one. Rank1 must use a scalar ray law, and rank0 a point atom; the proof handles both, while the frozen two-case compiler correctly operates only on its verified rank2 input.

The rational sector Laplace formula has the correct determinant, factorials and powers. It is first justified at sufficiently large positive transform parameters, where all signed terms are integrable. Summing and canceling the whole rational expression then comparing with the original product gives a valid exact identity. Unit mass at the origin follows after whole-expression cancellation and from the probability construction. No positivity is inferred from partial fractions or from finitely many samples.

The m2 triangular special identity and analytic transform formulas in the proof are consistent, including covariance `1/3600`, the lambda0/2 removable limits, and the stated triangular truncation tail. No kernel or special-function value was evaluated here.

## Generator preflight

Source lines47–57 implement exactly the determinant-bound interval selection above. Lines77–86 preserve actual columns, select a nonzero pivot, orient the sector and sort normalized directions counterclockwise. Lines94–102 implement both eta branches and merge only equal exponent covectors. Lines58–74 verify the whole exact rational Laplace identity and normalization after every stage. There are exactly seven intended checks: two stages for m3 and five for m6. The release pins, one-process/one-thread setting, nice minimum, no automatic restart and exclusive output creation are appropriate.

The resource caps are safe even before running the frozen cases. For m3, one old sector and one old term feed at most two new sectors, giving at most four endpoint contributions/final terms. For m6, the fan has one sector throughout: each old term yields at most two endpoint terms, so the stage term counts are bounded by 1,2,4,8,16 and the endpoint total by 30. Thus at most 34 literal endpoint contributions occur across both cases, far below100000; no stage can approach10000 terms. Degrees are at most4 and sectors at most2 for these cases.

Non-blocking bookkeeping limitation, separate from the mathematical PASS: line97 counts one resonant endpoint difference, which undercounts literal raw endpoints by one for that contribution, whereas the nonresonant branch counts its two endpoints separately. Treat the saved counter as generated combined expressions, or annotate this convention when reporting it. The independent bound34 controls literal endpoints regardless, so this does not invalidate the frozen release cap. It should be standardized before using this counter to size a future large batch.

## Actual independent verification

`registration.json` was written before the exact checks. `check.py` ran once through its symbolic jobs with SymPy1.14.0, nice10, one-thread environment; wall0.660054459 seconds, CPU0.394527 seconds. The exact m2 Jacobian/density, mass1, means1/20 and1/60, mixed moment1/900, covariance1/3600, and joint Laplace `1800/((s+30)(s+t+60))` all passed. Generic resonant/nonresonant antiderivatives were checked through degree4; the endpoint-comparison annihilation identity was checked symbolically. This is exact symbolic evidence supporting the source audit, not a numerical bound.

An earlier read-only preflight stopped before importing SymPy because this reviewer accidentally omitted the leading250 from the supplied manifest hash. The transcription was corrected before the single symbolic execution; the stop is recorded in `preflight_failure_receipt.json`.

Actual new m3/m6 density generations0; kernel/special-function values0; native calls0; old6D replays0. At completion none of `next_parent_release.json`, `next_started.json`, `next_result.json` or `next_exit_receipt.json` existed in the author packet. Main retains the release decision and later result audit. Main should carry the independent literal endpoint budget <=34 and this counter interpretation into the release and result report; do not label the program counter a literal raw endpoint count.

All writes belong only to `audit/upper147_review` in the independent no-hardlink clone. The source and main clones were read only. No commit or push was made.
