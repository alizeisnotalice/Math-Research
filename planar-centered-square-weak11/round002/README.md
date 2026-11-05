# R2 public evidence package

## Final project mapping

After integration, place the report at `proofs/round002-report.md`, the unchanged audit at `reviews/round002-audit.md`, this package's `verification/` at `round002/verification/`, `strict-interval.json` at `round002/strict-interval.json`, and `rounds/002.md` at `rounds/002.md`.

The proof is a derived public copy of frozen source SHA `5d9670c40f9630e0f873d0e3fb5a974a8a5b5260eb8f220b7031ea58bc07526e`: only status/path wording and references to non-public, unverified exploratory material changed; formulas and certified claims did not. The review (SHA `394d17e3c677f5b5132a8a6cc7ea4b3ab949819bd73cfdffcbd69233cc7f08a3`) and raw candidate input are byte-for-byte copies. Source/stage hashes and a summary of the non-mathematical transformations are in `verification/source-bindings.json`.

## Project-root replay commands

R1 already publishes both independent area checkers. This package does not copy or modify them. From the integrated project root:

```sh
python3 -I checks/lower/exact_check.py round002/verification/candidate002.json --Dr 1099511627776
python3 checks/lower/lower_check.py round002/verification/candidate002.json 1099511627776 4628923691329725842880531010053812883025951048952102577577/2718525473750034428308996140604808231510448537600000000000
python3 round002/verification/candidate002-transfer.py && wolframscript -file round002/verification/candidate002-scalar.wl
```

The second command uses the exact historical claim recorded in `candidate002-run.json`. That file and `.check2.json` preserve original machine paths as provenance; use the project-root-relative commands above for replay. The original area runs and their actual PIDs, exit codes, argv, and timings are preserved in the two core logs and run record; no area batch was repeated for this package.

## Scope

The two exact area implementations agree on the candidate's rational periodic lower certificate; the independent audit signs the general periodic-witness to finite strict-L1 transfer. The result is \(C_2\ge R_2>1.7027531463\), while R1's unchanged upper proof gives \(C_2<3.590867754339\). `strict-interval.json` records the exact endpoint gap and ratio. No optimality, exact \(C_2\), or attainment claim follows.

Publication relocation was additionally checked from a different working directory using the Fraction and Mathematica scalar scripts; the new actual run receipt is `verification/publication-scalar-replay.json`. Neither exact area batch was repeated.
