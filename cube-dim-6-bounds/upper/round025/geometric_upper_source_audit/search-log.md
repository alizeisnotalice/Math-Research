# C025 bounded source audit, 26 October 2026 (Asia/Shanghai)

Question: does an explicit strict improvement of the six-dimensional,
fixed-Lebesgue, centered axis-parallel cube constant follow from an
available covering refinement or from Lin's planar construction?

Scope fixed before search: four new queries, three primary documents,
20 minutes, eight named exact geometric checks. English web search used
web.run; no scholarly database or citation index was queried. Search hit
counts were not supplied and are not inferred. No search fragment is
used as theorem evidence. The four queries below were all executed on
2026-10-26 local date (2026-10-25 UTC):

1. `centered cube maximal weak type upper bound strict improvement 2^d Lebesgue`
2. `Besicovitch cubes strict gap covering maximal Lebesgue upper constant`
3. `Yongxi Lin "Two partial balayage principles" weak-type estimates covering`
4. `centered Hardy Littlewood cubes upper bound Stein Stromberg explicit constant`

Searches located the three author documents below and secondary discussions
of dimension-dependent estimates. The latter were not admitted as theorem
sources. The four-query search did not produce a verified explicit bound
below 64 for this problem; it does not establish that none exists.

The named historical `literature/search-log.md` was absent in the main
repository search. Existing README.md, primary_source_receipt.json and
round005 receipts were read instead. Their exact paths/hashes were locked
in input_manifest.json before this batch. Historical Aldaz Corollary 3.5
and Lin's planar scope are carryover background, not new progress.

## Admitted current primary documents (three, not three separate papers)

**D01** — Yongxi Lin, *Centred maximal constants for squares and Euclidean
balls*, author-maintained repository README.
URL: https://raw.githubusercontent.com/CoolRmal/centered-maximal-constant/main/README.md
Actual reading: title/authorship, lines 1–44, 108–131, 155–194, 252–257.
The relevant theorem table is lines 118–126 and the planar-kernel mechanism
is lines 157–174. Read scope: fulltext, partial sections only; not the
whole repository, no fresh Lean build. Hash equals the historical Lin
cache exactly. Identity/title/author/content verified; selected-version
publication/release date is not stated. Therefore bibliographic status
`pending_verification` for the year, no new BibTeX entry.

**D02** — Yongxi Lin, *The upper bound c_d <= 2^d*, file
`CenteredMaximal/UpperBound.lean` in the same repository.
URL: https://raw.githubusercontent.com/CoolRmal/centered-maximal-constant/main/CenteredMaximal/UpperBound.lean
Actual reading: lines 1–20, 47–59, 61–133, including the disjoint-selection,
center-cover, volume and mass summation, and limit arguments. Scope:
fulltext of the listed source portions, not dependency verification.
Copyright header gives 2026 and author; this is not silently equated with
a publication/release date. Bibliographic year remains pending; no new
BibTeX entry. The formal statement stays at 2^d. No execution of Lean.

**D03** — Yongxi Lin, *Two partial balayage principles for weak-type
estimates*, author's web article.
URL: https://coolrmal.github.io/articles/two-partial-balayage-principles/
Actual reading: title/author HTML lines 6–20 and 60–65; abstract/table
74–112; Sections 3.1–3.2 lines 366–455; Section 4.3 lines 525–555;
reference identities lines 723–744. Scope: fulltext, these sections only.
Principle 2 assumes an appropriate decomposition and admissible positive
kernel; its proof is explicitly an outline. Section 4.3 concerns planar
squares. We did not independently verify all analytic steps, fetch its
linked manuscript, or promote the additional AI-proposed numerical bound.
No publication date is supplied; HTTP Last-Modified and footer copyright
are not publication dates. Bibliographic status `pending_verification`
for the year, no new BibTeX entry.

D01 and D02 are two documents within one project, not distinct independent
confirmations. D03 links that project but is a separate author article.
This prevents double-counting corroboration. Exact bytes, URLs, retrieval
times, HTTP metadata and SHA-256 are in source_receipt.json. Mutable URLs
are pinned by those hashes locally, not claimed immutable upstream.

## Evidence and exclusions

- D01's formalized cube table has general 2^d and planar 3.879; its stronger
  higher-dimensional entries concern Euclidean balls. No six-dimensional
  cube theorem is inferred from the ball entries.
- D02 uses disjoint witness cubes and (1+tau)-enlarged cubes to cover their
  centers. Union-volume and mass inequalities are where the audited proof
  loses information. The new loss identity in theorem.md is our derivation.
- D03's maximal principle remains an admissible-kernel route; it supplies
  no ready six-dimensional covering refinement in the sections read.
- The retained historical Aldaz receipt distinguishes the supremum over
  measures from one fixed Lebesgue measure. It does not exclude a strict
  fixed-Lebesgue improvement. It is not counted as a fourth fresh source.
- No claim is made about the complete literature, publication priority,
  source correctness outside the read portions, or retraction status.

## Bibliography and stopping

references.bib carries one previously verified, directly relevant Aldaz
arXiv v2 entry, explicitly marked historical. The three current documents
are linked and hashed above, but missing publication dates prevent their
export under this skill's verification rule. No date was invented merely
to populate BibTeX. No fourth source or fifth query was opened to fill it.

The source audit is complete within its cap. The mathematical result is a
proved local/support-conditioned gap and a refuted overlap lemma, not a
new C6 bound. The next acceptance object is the explicit foreign-mass
selection lemma (L6) in theorem.md; no additional computation is launched.
