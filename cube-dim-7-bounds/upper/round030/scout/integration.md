# R030 numeric integration scout (read-only)

Date: 2026-10-27. Reviewed only the main checkout’s `upper/round030/compute` and `audit/round030/numeric` files. No scripts were run, no producer files were changed, and no new theorem search or MATLAB work was done.

## N7 result and smallest reviewable file set

The reported N30B strict threshold is

\[
U=12277313907857/500000000000=24.554627815714,
\]

with \(n=7\), \(t=3391749/524288\), and \(\sigma=2081927/6783498\). `compute/N30B.md:1-3,10,15` says the new upper bound follows **assuming** the general-\(\sigma\) two-region theorem in `geometry_theory/REPORT.md` §3.1.3; the exact certificates improve on the old same-\(t\) bound and the full round028 Arb lower endpoint. This scout did not inspect that theorem report, so it remains a required mathematical premise for final integration.

Suggested minimal whitelist for a reviewer to check the new N7 numerical bound without preserving search and exploratory diagnostics:

1. **Theorem and claim:** `cube-dim-7-bounds/upper/round030/compute/N30B.md` plus `cube-dim-7-bounds/geometry_theory/REPORT.md` §3.1.3 (the latter is expressly referenced by N30B but was outside this read scope).
2. **Producer arithmetic:** `upper/round030/compute/n30b.py`, which implements rational/log and integer-square-root enclosures and imports `revision_v3/numeric/verify.py` (`log_bounds`, `floorq`, `ceilq`). To reproduce its exact fixed-parameter certificate, preserve the `n=7` row of `n30b_certificates.json`; the full JSON also records checks and source hashes.
3. **Independent exact audit:** `audit/round030/numeric/independent_n7.py` and `independent_n7.json`, plus its exact comparison input `upper/round028/cap_candidate/result.json`. The auditor recomputes N7 rational intervals with a fixed 500-term log series and standard-library `Fraction`/`isqrt`; it does not import producer code. `audit/round030/numeric/PORTABILITY.md` documents the portability edit.
4. **If retaining the finite-a realization:** `upper/round030/compute/n30b_finite.py` and `n30b_finite_n7.json`; the independent follow-up is `audit/round030/numeric/finite_n7.py` and `finite_n7.json`. The audit script consumes `independent_n7.json` and checks the finite guard/cost using rational arithmetic.
5. **Shared scalar helper for producer reruns:** `revision_v3/numeric/verify.py`. Full regeneration of the multi-dimension certificate additionally uses `n30b_search.py`/`n30b_search.json`, `n30b_certify.py`, and `revision_v3/numeric/table.json`; these are not needed to audit the frozen N7 row or rerun the independent N7 audit. Search/stationarity traces, G30A root-isolation artifacts, the withdrawn L30A diagnostic, and other-dimension finite results are not needed for the N7 certificate.

## Runtime and portability

- Producer documentation records Python 3.15, SciPy 1.18.1, and mpmath 1.3.0 (`compute/README.md:3`). For the N30B exact arithmetic path the substantive operations use Python standard-library `Fraction`, `math.comb`, `math.isqrt`, `json`, and `pathlib`; `n30b.py` imports mpmath at module load, and `n30b_certify.py` uses mpmath for its additional positive-integral diagnostic. The N7 independent scripts use only the standard library. SciPy is not imported by the reviewed N30B N7 scripts and is not needed to check their rational certificate.
- `N30B.md:21-27` hard-codes `/Users/zhengzhihao/Desktop/cube-finite-dim-scratch/venv/bin/python`; this makes the printed commands host-specific. `compute/README.md:12-13` has the same absolute interpreter path and directs users to `stationarity.py`, which is absent from the reviewed `compute` directory. The certificate result itself uses project-relative inputs; the independent auditor’s file-location changes are described as project-relative in `audit/round030/numeric/PORTABILITY.md:1`. A portable reproduction should invoke the scripts with an available Python 3 interpreter and required mpmath, from the repository root (and install SciPy only if reproducing unrelated diagnostics).

## Contradictions and withdrawn material

- **Live stale parameter-mismatch text remains in producer source.** `N30B.md:17` explicitly retracts the old claim: \(\sigma=5/(3t)\) gives \(t\sigma=5/3\), and (R=7/4) is a valid separate guard because (R>t\sigma). However, `compute/n30b.py:37-45` still says at lines 39-40 that this requested pairing is inconsistent, and line 49 emits `parameter_mismatch: "sigma=5/(3t) implies R=5/3, not 7/4"`. The current `n30b_frozen.json:3` and `n30b_frozen_guard.json:12` correctly mark the old mismatch as retracted, and `n30b_frozen_guard.py:2,9-13` implements the valid guard. If `n30b.py` is rerun, it would overwrite `n30b_frozen.json` with its stale description; do not use that frozen routine/output to support or exclude N30B.
- **L30A exclusion is properly withdrawn.** `compute/L30A.md:3` says its old exclusion is invalid due to the wrong identity and coefficient; it retains the old script/JSON only as a record and directs use of `l30a_corrected.py`/`l30a_corrected_result.json` for the corrected scalar diagnostic. Neither old L30A exclusions nor its exploratory diagnostics belong in the N7 whitelist.
- **Ancillary README command is stale.** `compute/README.md:12` names missing `stationarity.py`. This does not invalidate the N30B certificate, but should not be presented as a reproducible command for the current directory.

No other N30B parameter-description inconsistency was found in the reviewed files. No conclusion here relies on the invalid L30A exclusion.
