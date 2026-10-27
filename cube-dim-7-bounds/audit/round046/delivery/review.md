# R046 delivery review

Read-only review of the compute report, registration, case summary, monitor, and symbolic check. No numerical experiment was rerun. No producer or main-clone file was modified.

## Identity and execution evidence

- Report: `/Users/zhengzhihao/Desktop/T/cube-r030-compute/cube-dim-7-bounds/upper/round046/compute/R046_REPORT.md` (SHA-256 `eb83cf8f00d7ffa211dd41e70537ffa5d2ccac6f146e6c3ef57726323a4f066c`).
- Registration: `registration.json` SHA-256 `6066e296172a1c7c9ea4b64d58c56669552ed482dd58f87fe64b4f0096f9fda5`; it pins `r046_row.py` hash `d6bed3752b35b17bdacad572ec3e1b4ed02ece2d23f669d28facc60402a2310f`, 15 workers, 80 decimal digits, 3,840 combinations, and a stop-on-failed-comparison rule. The report and registration timestamps follow the source file timestamp; the report was generated later. `symbolic_check.py` is a separate post-batch exact algebra check with its own script hash, not part of the registered numerical worker source.
- `summary.json` agrees with registration: 15 distinct task PIDs, each with 97 direct, 97 mean, and 97 second-moment integral calls, totaling 4,365 calls for 3,840 `(h,r)` combinations (identical `r` values are cached). Summed task CPU is 36.565743 s; batch wall time is 2.760134 s, or about 13.25 CPU cores averaged over the elapsed batch interval.
- The single `monitor.jsonl` sample at 0.098 s records 15 workers and 0.0% aggregate instantaneous CPU. The tasks ended before the later five-second sample. The integral CPU/wall ratio does not establish an instantaneous full-core snapshot or sustained saturation.
- Identity hashes: `summary.json` `0760d918e1876fb0341c3d35f8777182760e669ea37d721976ce02cd4f70df44`; `monitor.jsonl` `3d378d5c1a2dcc553c36e8286cdc2d347580af5827e442e4eaf564901566eafd`; `symbolic_check.json` `62dbe3e0217efad962df141128857e6f9f34aa1b4cfa524eec967ee3e3d8a2c5`; numerical source `r046_row.py` `d6bed3752b35b17bdacad572ec3e1b4ed02ece2d23f669d28facc60402a2310f`; symbolic source `symbolic_check.py` `f7fdbece7f1394c7cecbaa595ca6ccb8e569fef564a621d9e2a267630dfc2940`.

## Mathematical and numerical scope

- The computation uses mpmath at 80 decimal digits to compare the defining piecewise integral with the proposed closed form, then checks finite-case face-mean and second-moment bounds. The maximum observed absolute direct/closed discrepancy is approximately `1.20e-81`; it is a finite-precision discrepancy, **not** a rigorous interval error bound.
- The finite grid's largest normalized mean deficit and second moment are approximately 1.967303 and 2.211449. These observations do not prove the proposed bounds 6 and 31 uniformly, sharply, or over all parameters.
- `symbolic_check.json` reports PASS for the identity defining `E(u)` and the algebraic face-mean cancellation at `r=R` under the stated normalized-domain assumptions. This does not establish the general uniform second-moment bound.
- The report describes the work as one-dimensional numerical diagnostics, not an interval proof or `C_n` theorem. No covariance grid or extra asymptotic experiment is reported. The run does not establish a new `C_n` bound or a seven-dimensional lower bound.

## Reproduction caveat

The report's commands use the absolute interpreter path `/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python`, outside the R046 source directory. That interpreter exists on this host, but the registration pins source and parameters, not the virtual-environment contents or package versions. Reproduction elsewhere therefore depends on that external venv or a separately provisioned compatible environment. This is an environment portability caveat, not evidence that the reported run did not execute.

## Conclusion

Counts, CPU/wall ratio, finite precision, and the report's stated limitations match the registration, summaries, monitor, and symbolic result. Do not describe the `1.20e-81` discrepancy as a certified error enclosure, the sampled maxima as global or sharp constants, the CPU ratio as instantaneous sustained full load, or this run as a new `C_n` result. No report correction is indicated from this review.
