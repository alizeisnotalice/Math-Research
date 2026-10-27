# R052 archival whitelist

Copy these files under `upper/round052/compute/`:

- `R052_REPORT.md`, `R052_REPORT_producer_original.md`, `ACTIVATION.md`, `ARCHIVE.md`;
- original frozen `registration.json`, `results.json`, `r052_integer.py`, `equivalence.json`, `verify_top64.py`, `verify_top64.json`, `monitor.jsonl`;
- `portable_check.py`, `portable_smoke.py`, `portable_smoke_result.json`, and `archive_checksums.sha256`.

`results.json` is approximately 336 KiB and is the sole exact candidate bundle; there is no compressed duplicate. It includes all 64 global and 15×32 worker certificates. The original per-chunk top lists and an independent chunk execution ledger were not retained. The occupancy argument establishes top-64 completeness only conditional on the frozen producer's per-chunk execution and retention behavior; retained checks and exact counts do not independently prove that every declared chunk ran exactly once.

The original producer report is retained byte-for-byte as `R052_REPORT_producer_original.md` because the independent audit records its original hash. `R052_REPORT.md` changes only the confusing wording around the retained-candidate verification command; the historical conditional theorem status is otherwise unchanged. Later accepted R051/R053/R055 results are described separately in `ACTIVATION.md` with proof and audit hashes; the original R052 run is not rewritten as having known them.

Also archive the separate independent retained-point audit under `audit/round052/`, and the accepted R051/R053/R055 proof and independent audit sources referenced in `ACTIVATION.md`. The producer's own Fraction check is not an independent author audit. The R052 auditor's scope is retained points and aggregate counts, with exhaustive execution still conditional.

Run `python3 upper/round052/compute/portable_check.py` in the copied project layout. It uses only standard-library Python and relative files, checks exact retained certificates and list-count arithmetic, and does not rerun the 7,372,800 tuples or modify original receipts. `portable_smoke.py` tests this minimal set in a clean temporary project tree; it is an archive-adaptation check, not a new numerical experiment.
