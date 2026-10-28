# R070 archive

Copy this whole folder, including manifest and temporary verification receipt.

- `finite_width.py`
- `finite_width_e.py`
- `h2.py`
- `original_chunks.json.gz`
- `pilot_monitor.jsonl`
- `pilot_summary.json`
- `production_monitor.jsonl`
- `production_summary.json`
- `registration.json`
- `verify_summary.py`
- `width_e_monitor.jsonl`
- `width_e_registration.json`
- `width_e_summary.json`
- `width_monitor.jsonl`
- `width_registration.json`
- `width_summary.json`

Pilot and production H2 use 160 decimal closed form and 256 decimal independent checks, with one-sided kink derivatives. width is explicitly a surrogate n=2^L, not theoretical L=ln n; width_e is the separately registered corrected n=exp(L) diagnostic. Width_e carries rho*(J−1) separately and centered variance, though tiny inverse corrections may still round at 90 digits. Both width stages use 32/64 point split quadrature; neither is a strict proof. Original chunk UTF-8 files are losslessly retained in one gzip. Cantelli estimates apply only to stated model.

Run python3 verify_summary.py to verify retained bytes/counts; it does not rerun integrations.
