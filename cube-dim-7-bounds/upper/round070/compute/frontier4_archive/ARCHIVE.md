# R070 conditional frontier archive

Copy this directory including manifest and temp verification receipt.

- `frontier4.py`
- `frontier4_monitor.jsonl`
- `frontier4_recheck.json`
- `frontier4_recheck.py`
- `frontier4_registration.json`
- `frontier4_summary.json`
- `original_chunks.json.gz`
- `verify_summary.py`

All 1280 original UTF-8 receipts are losslessly stored in one gzip. Float screening is not a proof; 120-decimal S64 candidate recheck is post-run, and other finite guards/operator hypotheses were not established by this batch. No continuous optimum or C_n bound follows. Run python3 verify_summary.py to check retained hashes/counts/top64 only.
