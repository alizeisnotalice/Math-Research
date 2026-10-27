# R057 archive

Copy this entire archive directory preserving relative paths, including `archive_manifest.json` and this file.

- `compact_receipts.json.gz`
- `corrected_v2/REPORT.md`
- `corrected_v2/monitor_pilot.jsonl`
- `corrected_v2/monitor_run.jsonl`
- `corrected_v2/r057_kinks_v2.py`
- `corrected_v2/registration.json`
- `corrected_v2/summary_pilot.json`
- `corrected_v2/summary_run.json`
- `monitor.jsonl`
- `postrun_precision_probe.json`
- `postrun_precision_probe.py`
- `r057_kinks.py`
- `registration.json`
- `summary.json`
- `verify_summary.py`

Initial 240 failures have generic assertions without line numbers: their individual causes are not established. Corrected 240 complete receipts are retained losslessly. The post-run one-point 160/220/256-digit precision probe is independently recorded after the original run; it explains one representative failure but cannot classify all 240. The earlier single-point original tool output was not retained as a raw receipt.

Run `python3 verify_summary.py` for retained hashes, counts and probe receipt; no numerical experiment is rerun.
