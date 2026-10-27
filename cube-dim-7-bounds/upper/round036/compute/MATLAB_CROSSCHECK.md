# R036 MATLAB exponent correction cross-check

The older T35C explanation for the 51200 L^-5/4 cap bound was incorrect: direct e^-u≤u^-2 yields 51200 L^-3/4. The revised analytical proof for z=log L≥128 is described in corrected `upper/round035/compute/T35C.md`. The underlying old Arb boxes were not rerun or modified.

The MATLAB script `r036_cap_crosscheck.m` was executed through the connected MATLAB tool under the coordinator's exclusive window. It used only task-prefixed variables, no clear/resetpath/path change or cwd change. Its exact symbolic polynomial values and 40-digit registered-point logarithmic comparisons are recorded verbatim in `matlab_result.txt`. This check is an independent computational diagnostic and does not replace the uniform analytic proof or interval arithmetic. The exclusive window was released after execution.

SHA-256: script `9fe43dba6532b5b0e9fde2e7f6f68a1972cad4e3a12bccd5cdae18fa81fe41c1`; recorded output `7b86537b63183674370d5b60d8d7db464cd2b6b22d56d6b1797e61656b1906d6`.
