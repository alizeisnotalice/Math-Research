# R044 portable archive whitelist

Copy the following files under `upper/round044/compute/`:

- `R044_REPORT.md`, `ARCHIVE.md`, `r044_compact.json.gz`, `archive_manifest.json`, `verify_archive.py`, `archive_verify_receipt.json`;
- `registration.json`, `summary.json`, `triple_retry_registration.json`, `triple_retry_summary.json`, `B_registration.json`, `B_summary.json`, `C_registration.json`, `C_summary.json`;
- `resource_monitor.jsonl`, `retry_resource_monitor.jsonl`, `B_monitor.jsonl`, `C_monitor.jsonl` (the last may be empty because the batch finished before its first sample);
- `r044_batch.py`, `r044_triple_retry.py`, `r044_b.py`, `r044_c.py`, `make_archive_bundle.py`, plus `archive_checksums.sha256`.

Also archive the R038 `phasea_compact.json.gz` under `upper/round038/compute/`. The portable checker resolves it by this relative project layout; it needs only Python standard library and does not read any private absolute path, original raw JSON or solver. R038's independent audit is the provenance for the cached geometry and clique lists; this producer-data replay is not a new independent reconstruction of those lists.

The R044 compact bundle preserves the original exact strings for all 76 primal/dual vectors, starting and global feasible q, first-32 best vectors and statuses, both B candidates, all 38 C vectors and all explicit cases without vectors. It maps each case to the original raw files' SHA-256 hashes. The 2432 and 6702 per-proposal logs remain in the private clone and are deliberately not archived; neither aggregate non-improvement nor floating MILP status proves optimality. The 76 original full per-case files are also retained privately. Original solver failures and retry registration remain represented in the summaries and report.

From the copied project layout run:

```
python3 upper/round044/compute/verify_archive.py
```

The checker validates all 76 rational primal/dual inequalities, payments and global-load comparisons, all saved triple and C vectors, against the exact audited R038 cache; it writes a receipt beside itself. It does not re-solve any LP, claim minimal support, or independently regenerate cube geometry.
