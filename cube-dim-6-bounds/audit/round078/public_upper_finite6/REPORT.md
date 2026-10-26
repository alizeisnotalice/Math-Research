# Round 078 — fixed six-dimensional public upper bound scan

Search date: 2026-10-26 (UTC; Asia/Shanghai workspace date also 2026-10-26).

## Result

This directed scan did not establish a new explicit numerical upper bound below the current public (C_6\le64) for the weak-(1,1) constant of the Lebesgue centered, axis-parallel, common-radius cube maximal operator on general (L^1(\mathbb R^6)).

The main lead is Stein–Strömberg’s original covering theorem. It applies to this cube family: for every bounded, open, convex, symmetric body (B\subset\mathbb R^n), the centered maximal operator over its dilates obeys a weak-(1,1) estimate of order (c n\log n), with (c) independent of both (B) and (n). The paper’s Theorem 1 states existence of an absolute (c), but does not state a numerical value for it. Its proof also leaves the covering constants unquantified at the point needed to produce a reproducible (n=6) value. Thus (6c\log6) cannot be compared numerically with 64 from the published statement, and no improvement is claimed. See the original paper’s Theorem 1 (print p. 260; proof, print pp. 261–264) and the detailed source receipt.

The same paper’s Theorem 3 improves the order to (c n) only when (B) is the Euclidean unit ball, using a heat-diffusion semigroup. This is a different averaging body from an axis-parallel cube, and its coefficient is again not given as a numerical constant. It cannot be substituted for the cube operator.

The current literature ledger’s public upper bound remains (2^6=64), from Aldaz, *Kissing numbers and the centered maximal operator*, arXiv:1811.10478v2, Theorem 1.2 and Corollary 3.5: the uniform centered constant for the ℓ∞ norm is (2^d), which applies to Lebesgue measure and centered axis-parallel cubes. The current receipt records its PDF URL, SHA-256 `d42a1d6518298227e2af47876761d47225e38b28e2f7726f91c0cc2498ab12e7` and reading range (PDF p. 2; pp. 10–11). This supplemental task did not download or re-read it and did not redo any Aldaz root calculation.

No explicit upper recurrence specialized to centered cubes in dimension was located in this bounded search. A cited monotonicity statement (c_{d+1}\ge c_d) is a lower-direction relation, not an upper recurrence. This report records a search limitation, not a claim that no better result exists.

## Scope and distinctions

- **Target:** (M_6 f(x)=\sup_{r>0}|Q(x,r)|^{-1}\int_{Q(x,r)}|f(y)|\,dy), with (Q(x,r)=x+[-r,r]^6), Lebesgue measure, every scale and general (L^1(\mathbb R^6)) input.
- **Applicable but non-numerical:** Stein–Strömberg Theorem 1 covers bounded open convex symmetric bodies, hence open axis-parallel cubes, but only provides (c n\log n) with unspecified universal (c).
- **Different body:** their Theorem 3 concerns Euclidean balls, not cubes; its (c n) coefficient is also unspecified.
- **Different objects excluded:** Besicovitch constants for arbitrary measures are used only through the centered-ℓ∞ result already in the current ledger; bounds for independent rectangles, strong (L^p), Euclidean balls, or lower bounds are not relabeled as an improved (C_6) upper bound.
- **Finite-dimensional recurrence:** no upper recurrence was established by the eight-query budget. The known monotonicity (c_{d+1}\ge c_d) points the other way.

## Original source actually read

Stein, Elias M.; Strömberg, Jan-Olov, “Behavior of maximal functions in \(\mathbb R^n\) for large \(n\),” *Arkiv för Matematik* 21, 259–269 (1983), DOI [10.1007/BF02384314](https://doi.org/10.1007/BF02384314). Publisher metadata gives publication in December 1983 and authors/volume/pages. A public full-text copy was downloaded from the Tsinghua archive and retained here. `sources.json` records both URLs, version, byte count and SHA-256.

Actual reading: the 11-page PDF text was inspected for the title/introduction and object definition (PDF pp. 1–2); Theorem 1 and its covering selection/constant estimates (PDF pp. 2–6, print pp. 260–264); Theorem 3’s separate Euclidean-ball semigroup argument (PDF pp. 6–7, print pp. 264–266); and the appendix’s subject matter and reference list (PDF pp. 8–11, skimmed only to distinguish the (L^p)/spherical material from the weak-(1,1) cube question). Relevant extracted positions are in the source receipt. The paper was accessed through the archived PDF and publisher metadata page, not via a secondary summary.

## Search log and stopping limit

Eight directed web queries were used; the query cap is exhausted. Queries and source-selection notes are in `search_log.md`. The source scan was bounded and targeted. No broad or exhaustive literature search was attempted. No conclusion of nonexistence is made.

No scientific kernel, numerical search, old certificate script or matrix/scalar re-computation was run. No MATLAB was used. The existing Aldaz (n=6) root certificate was not rerun. Only this directory was written; no staging, commit or push occurred.
