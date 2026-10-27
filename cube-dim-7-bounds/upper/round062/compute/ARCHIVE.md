# R062 compact archive

Copy this entire directory preserving names, including `ARCHIVE.md`, `archive_manifest.json`, and the post-run `temp_verify_receipt.json`. Files:

- `candidate_cert.py`
- `candidate_certificate.json`
- `candidate_gate.json`
- `candidate_gate.py`
- `frontier.py`
- `frontier_monitor.jsonl`
- `frontier_registration.json`
- `frontier_summary.json`
- `frontier_v2.py`
- `frontier_v2_monitor.jsonl`
- `frontier_v2_postrun_manifest.json`
- `frontier_v2_registration.json`
- `frontier_v2_summary.json`
- `frontier_v2_top_cert.py`
- `frontier_v2_top_certificate.json`
- `original_chunks.json.gz`
- `verify_summary.py`

The two sets of original chunk UTF-8 bytes are retained losslessly in one gzip container; source registration, original summaries, exact fixed 128-layer candidate, post-run gate, resource monitors and post-run top certificate are retained. The post-run source-to-chunk mapping is not a field of the original worker receipts.

Pilot: 4096 tuples / 245760 layers; v2: 131072 tuples / 50331648 layers, 0 overlap by exact Fraction registration. V2 wall 136.02 s, summed worker CPU 1846.22 s (average 13.57 cores), sampled peak 1458.8% CPU and 458416 KiB RSS. The numerical top is a finite model candidate, not a continuous optimum or C_n theorem.

Run `python3 verify_summary.py` to check archive bytes, original chunk hashes, counts, retained top lists and exact recorded arithmetic. The check does not recalculate the omitted layer terms, independently verify original computational coverage, or re-prove the operator theorem.
