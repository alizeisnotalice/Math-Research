# R067 compact archive

Copy the entire directory, including this file, archive_manifest.json and temp_verify_receipt.json.

- `candidate.py`
- `candidate_certificate.json`
- `monitor.jsonl`
- `original_chunks.json.gz`
- `phi_diagnostic.py`
- `ratio_probe.json`
- `ratio_probe.py`
- `registration.json`
- `summary.json`
- `verify_summary.py`

256 original chunk UTF-8 files are retained byte-for-byte in one gzip. Producer registration/source/summary/monitor and nine-gate scalar source/certificate retain original identities. ratio_probe files are a separate post-run 1536-state H/J diagnostic, never attributed to the production run. The 512 baseline cases each include a mean check and second-moment/variance computation; they are not 1024 separate cases.

The closed-form production used 160 decimal digits; independent checks used 256 decimal digits (>512 bits). Capture-onset derivatives were evaluated one-sided; this finite diagnostic does not prove an everywhere derivative or global theorem. Original wall 70.40 s, worker CPU 953.83 s, sampled peak 1462.2% CPU and 3765840 KiB RSS (~3.59 GiB).

Run `python3 verify_summary.py` to verify retained bytes/hashes/counts/recorded fields only. It does not rerun any integrations or derivatives.
