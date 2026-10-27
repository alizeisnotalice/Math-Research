# R047 portable archive whitelist

Copy under `upper/round047/compute/`:

- `R047_REPORT.md`, `ARCHIVE.md`, `archive_summary.json.gz`, `archive_manifest.json`, `verify_archive.py`, `archive_verify_receipt.json`, `archive_checksums.sha256`;
- `registration.json`, `summary.json`, `monitor.jsonl`, `cross_face.py`, `micro.json`, `check_moments.py`, `check_moments.json`, `make_archive.py`;
- coordinator-produced `r047_primary_matlab_aggregation.m` and `r047_primary_matlab_aggregation.json`, copied unchanged with hashes. They document the coordinator's actual exact finite-law pinned-face aggregation check for n=2,…,7, not a continuous bound or an n=7 maximal-operator experiment. The original MATLAB script contains a private absolute receipt path and is retained as a historical source, **not** as the portable verifier; do not execute it as part of archive verification.

The 15 full raw case files remain in the private clone and are not archived. The compact summary preserves each case's smallest and largest stored numerical observation, registered parameter triple, integral/resource counts and original raw SHA-256. The producer hash fixes the original code. `make_archive.py` is retained for provenance but needs the omitted raw files to reconstruct the summaries.

Run `python3 upper/round047/compute/verify_archive.py` from the copied project layout. The standard-library checker reads only archived paths and checks stored extrema, parameter membership, basic numeric bound arithmetic, recorded verification statuses and source hashes. It cannot independently establish maximality over omitted raw rows, rerun the cross-moment integrals, or turn numerical values into rigorous interval bounds. The local reproduction command in the research report is not required by the portable verifier.
