# R058 archive

Copy this entire archive directory preserving relative paths, including `archive_manifest.json` and this file.

- `certify_model.py`
- `compact_receipts.json.gz`
- `model_certificate.json`
- `monitor.jsonl`
- `r058.py`
- `registration.json`
- `summary.json`
- `verify_summary.py`

128 original chunk receipts retained losslessly compressed, including top32; original SHA identities and summary top64 preserved. Omitted tuples were not stored; no independent full-grid re-evaluation or continuous optimum.

Run `python3 verify_summary.py`. This post-run checker verifies retained hashes, record counts and retained certificates only; it does not redo the underlying numerical experiments. The archive generator and checker are post-run artifacts, not independent authors or original run evidence.
