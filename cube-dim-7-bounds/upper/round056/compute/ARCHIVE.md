# R056 archive

Copy this entire archive directory preserving relative paths, including `archive_manifest.json` and this file.

- `R056_REPORT.md`
- `compact_receipts.json.gz`
- `monitor.jsonl`
- `r056_inverse.py`
- `registration.json`
- `summary.json`
- `symbolic_identities.json`
- `symbolic_identities.py`
- `verify_summary.py`

60 compact receipts retain original hashes, counts, PID and errors; 16 MB original per-state rows omitted, independent Luna read those originals.

Run `python3 verify_summary.py`. This post-run checker verifies retained hashes, record counts and retained certificates only; it does not redo the underlying numerical experiments. The archive generator and checker are post-run artifacts, not independent authors or original run evidence.
