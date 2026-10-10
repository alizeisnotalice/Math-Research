# 2026-10-11 current release snapshot

This directory contains the public, machine-readable metadata snapshot for the 64 revised math-research Skills in `../../skills/`. It is a metadata and provenance release, not a theorem certification. The local evidence bundle is distributed separately; this repository change adds no paper PDFs or full manual prompt/response records.

## What the snapshot says

- The current index covers 64 Skills. Its status is `FINAL_WITH_DOCUMENTED_LIMITATIONS`; the package validator performs structural checks and explicitly does not establish mathematical proof or successful Skill runtime behavior.
- The canonical claim catalog contains 444 ledger-derived rows. The rows preserve their source, assumptions, proof/evidence axes, review state, and known gaps. Catalog inclusion is not claim acceptance.
- The dependency graph contains 64 Skill nodes, 2,791 source-version nodes, and 3,089 declared-source edges. It records references in source manifests; it is not a theorem-dependency graph.
- The current main-source screen covers 93 exact-version sources across 64 Skills (5,952 source–Skill decisions): 5,736 unrelated at method scope, 192 adjacent, and 24 direct. There are zero claim transfers and zero whole-paper proof acceptances in this screen. These are bounded method-scope decisions, not theorem proofs.
- The additional-source snapshot covers 85 sources × 64 Skills (5,440 decisions) and keeps partial reads, pending decisions, version gaps, and unavailable source bytes explicit. Read-extent buckets are 60 full-page records, 1 historical-version full read with the requested version unavailable, 5 records without an exact-SHA read card, 17 records whose extent is not verified full, and 2 preview-only records. This is not 85-source full-text coverage.
- The legacy 18-Skill crosswalk has 148 claim-heading rows with candidate retrieval, but zero confirmed exact current-claim mappings. Its historic audit remains pending current-input binding; lexical candidates are not semantic equivalences.
- The fixed application matrix has 128 slots. Current prompt coverage is 128/128; 24 recent records are explicitly manual and unblinded, and 104 older records retain unknown/null external-call fields. No native Skill-task execution is claimed. Historical prompt gaps remain unrecovered.
- Five source records lack verified local PDF bytes. Their resolution is recorded as unavailable/record-only, not as a successful paper lookup.

Counts are separate by evidence type. Source fit, ledger status, examples, installed runtime, and proof verification must not be conflated. See the JSON records for exact bindings and scope.

## Files

- `index-64.md` and `index-64.json`: current Skill inventory and status.
- `claim-catalog-current.json`: canonical 444-row claim ledger derivation.
- `dependency-graph-current.json`: declared source links.
- `main93-current64-source-scope-report.json`: 93-source × 64-Skill method-scope report.
- `additional85-current64-source-scope-read-overlay.json`: supplementary source scope and read-extent report.
- `legacy18-current-crosswalk.json`: bounded legacy retrieval crosswalk and pending semantic mapping status.
- `files.json`: exact SHA-256 and byte count for this metadata directory (the manifest does not hash itself).
- `skills-files.json`: exact SHA-256 and byte count for each file in the public 64-Skill tree.

The complete local evidence package contains source files and audit material that are not mirrored to this public repository. Do not infer that the repo snapshot contains all papers or local caches.
