# R046 portable archive whitelist

Copy under `upper/round046/compute/`:

- `R046_REPORT.md`, `ARCHIVE.md`, `archive_summary.json.gz`, `archive_manifest.json`, `verify_archive.py`, `archive_verify_receipt.json`, `archive_checksums.sha256`;
- `registration.json`, `summary.json`, `monitor.jsonl`, `r046_row.py`, `symbolic_check.py`, `symbolic_check.json`, and `make_archive.py`.

The 15 full raw case files remain in the private clone and are not archived. The compact summary keeps the registered parameters, each case's worst observed rows, actual integral counts/resources and all raw SHA-256 hashes. The source code hash fixes the original producer; `make_archive.py` is retained for provenance but requires the unarchived raw files to regenerate the summary.

Run `python3 upper/round046/compute/verify_archive.py` from the copied project layout. This standard-library checker reads only archived files and checks the stored observations, parameter membership, numerical bound arithmetic and provenance hashes. It cannot independently prove that a saved observation is the maximum over omitted raw rows, nor does it rerun quadrature or provide rigorous interval enclosures. The reproduction command in the research report uses the local environment of the original production run; it is not required by the portable verifier.
