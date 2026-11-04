# Round 29 truncation audit

Run from this directory:

```sh
python3 round29_truncation.py fresh
python3 round29_truncation.py replay
wolframscript -file round29_truncation.wl
```

Run `round29_truncation.m` through MATLAB.  Python writes a temporary local cache only to demonstrate fresh/replay; it is intentionally not part of the archive.  `round29_summary.json` is the compact retained output.

The continuous Pareto density uses the exponent `q - 1/gamma`.  Both rational split diagnostics and the finite `q=1/gamma` logarithmic partial-moment check run on every Python invocation.
