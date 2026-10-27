# R065 archive

Copy this directory, including manifest, this file and temp_verify_receipt.json.

- `candidate.py`
- `candidate_certificate.json`
- `candidate_initial_sign_error.json`
- `candidate_initial_sign_error.py`
- `derivative.py`
- `derivative_monitor.jsonl`
- `derivative_registration.json`
- `derivative_summary.json`
- `original_chunks.json.gz`
- `verify_summary.py`

Original 128 worker chunk UTF-8 bytes are losslessly held in one gzip. Initial scalar source had a sign transcription error, falsely failing tau; both that version and the corrected certificate are preserved. The original error is not a mathematical counterexample. Independent numerical checks used 256 decimal digits (>512 bits), not 512 decimal digits.

Run `python3 verify_summary.py` to check retained data and scalar identities. It does not rerun integrals or establish global bounds.
