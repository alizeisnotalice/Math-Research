# Round 29 computation archive

This is the corrected v2 truncation audit for the unbounded cube angular profile.  The continuous angular Pareto density is \((n-1)(n\gamma)^{-1}a^{-1/\gamma-1}\), so its truncated \(q\)-moment has exponent \(q-1/\gamma\), not \(\gamma q-1\).  The earlier scratch v1 finite-truncation output using the wrong exponent was discarded and is not archived.

The retained archive is exactly `src/round29/{README.md,round29_truncation.py,round29_truncation.wl,round29_truncation.m,round29_summary.json}` plus this report.  The compact JSON omits the reproducible 96-row table, raw logs and caches.

Actual commands from `src/round29` were:

```sh
python3 -m py_compile round29_truncation.py
python3 round29_truncation.py fresh
python3 round29_truncation.py replay
wolframscript -file round29_truncation.wl
```

MATLAB R2026a statically checked and executed `round29_truncation.m`.

Python ran 24 parameter triples at 120 decimal digits, with 15 actual fresh worker PIDs 4686--4700 and zero failures; replay had 24/24 cache hits.  It executes two exact `Fraction` split diagnostics, \((n,\gamma,q,A)=(2,1/2,0,2)\) and \((2,1/2,1,4)\), on every run, and also checks the finite logarithmic partial branch \((2,1/2,2,2)=1/2+\log2\).  Wolfram 15.0 passed its moment checks, 216 actual \((n,\gamma,c;q,A)\) split checks, both rational diagnostics and the logarithmic branch.  MATLAB independently evaluates 96 finite-radius values and 216 Pareto density-split controls; its reported maximum quadrature split error was \(9.95\times10^{-14}\), below the \(10^{-9}\) assertion threshold.

These are finite-truncation implementation checks after volume normalization.  They do not assert a joint \((\tau,R)\) limit or replace the analytic rearrangement proof.
