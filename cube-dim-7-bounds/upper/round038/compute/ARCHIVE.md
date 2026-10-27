# R038 archival whitelist

Recommended `upper/round038/compute/` files:

- `R038_REPORT.md`, `ARCHIVE.md`, `archive_checksums.sha256`, `phasea_registration.json`, `phasea_summary.json`;
- `phasea_compact.json.gz`, `compact_manifest.json`, `replay_compact.py`, `replay_receipt.json`;
- `phasea.py`, `pair_balance.py`, `joint_core.py`, `joint_geometry.py`, `phasea_cert3.py`, `phasea_cert3_registration.json`, `phasea_cert3_results.json`;
- the three `lp_cert_phasea_*.json` files, which retain exact primal/dual vectors; their corresponding geometry, budgets, U, cliques and q-before live in the compact bundle and final q is reproducibly reconstructed;
- `make_compact.py` as a provenance record. It requires the 79 raw source JSONs still retained in the private clone and cannot regenerate them from the archived bundle.

Also retain the soft-agent R038 SPEC and the independent audit reports, code, and essential audit result summaries under their respective project paths. Preserve the R037 historical geometry generator and its provenance because `phasea.py` references it; replay itself is standalone and does not invoke that generator, original LP or old atomization. The 79 large raw result JSONs remain in the private clone but are intentionally not on this archive whitelist. Their SHA-256 values and final-q canonical hashes are in `compact_manifest.json`; the compressed bundle is not a byte-for-byte substitute for those raw files.

Verify without SciPy or original raw files:

```
python3 replay_compact.py
```

Use `python3 replay_compact.py --terminal` only if a fresh terminal nonedge scan is wanted; that optional scan is not needed for the current archived receipt and does not reproduce the original search ranking. The full independent audit already checked all 50 terminal conditions.
