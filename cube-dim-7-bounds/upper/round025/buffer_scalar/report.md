# R025 free-buffer finite scalar certificate

The best of the three registered optimized-threshold rows is **t=11/2**, giving the strict centered-cube weak-type bound

    C7 < 29059800754803 / 10^12 = 29.059800754803.

The underlying 512-bit interval is approximately 29.0598007548021301377, strictly below R024's full candidate interval near 46.6442895597215742010. This selects only among the three specified t values; it does not claim a global t optimum or a lower bound on the optimal constant.

The separately frozen theorem, independently audited, proves for every t>1 the original centered convex-body operator bound

    Cn ≤ [(1+1/t)^(n/2)+√(1+n log t)]².

It retains arbitrary L¹ input, the same centered body and every positive common radius. The finite computation substitutes n=7; it does not itself replace the selection proof or the passage to the limiting buffer parameter.

| t | optimized-threshold expression, approximate | strict 10^-12 upper grid |
|---:|---:|---:|
| 5 | 29.1084362054948297 | 29.108436205495 |
| 11/2 | 29.0598007548021301 | **29.059800754803** |
| 6 | 29.1080175070108745 | 29.108017507011 |

The fourth predeclared calculation uses rational parameters t=5, τ=20/3 and gives approximately 29.1086423516626165, strictly below 2233/75. It shares the log 5 evaluation with the first row. All four 512-bit outward dyadic intervals are in terminal.jsonl; each optimized row is strictly separated from the other two, and its ceiling is floor(U·10^12)+1 over 10^12.

The scalar run made exactly four evaluations, with three log and six square-root calls, and no parameter search. Wall time was 0.00130 s, CPU time 0.00130 s, and peak process RSS 30,392,320 bytes. The checker repeats the same four scalar inputs, verifies both frozen proof hashes, the R024 complete interval, the strict ordering and grid rounding. It does not establish the operator theorem or claim a continuous optimum. No extra t, old baseline, optimizer, MC, phase, MATLAB, commit or push was used.
