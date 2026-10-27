# R052: exact finite conditional-parameter search

The best candidate in this finite family is

\[
 a=73/8,\quad s=293/16,\quad d=1023/1048576,\quad
 c_{\rm raw}=\frac{6468305130094850733787980147}{370780099533215320770594052505600}
 \approx1.7445124854969205\times10^{-5}.
\]

This is a coefficient of the **conditional scalar model** provided to this batch. It is not by itself a bound for C_n, an improvement of a complete operator estimate, or a continuous-parameter optimum. The general R051 argument, its independent audit and comparison of all costs are required before any such conclusion. Each parameter must remain fixed when applying an eventual-n theorem.

The finite list consists of a=(32,…,256)/8, ρ∈{1,2,4,8,16,32,64,128}/16, s=2a+ρ, and d=k/2²⁰ for k=1,…,4096: exactly 7,372,800 tuples. The recorded run processed all 1,800 deterministic chunks once. It counted 3,975,632 valid tuples and 3,397,168 rejections by d<1/(32a); there were no σ or bracket guard rejections. It used no grid outside the registration and stores only retained candidates, counts, hashes and resources.

All ranking uses integer cross multiplication. In this registered domain M=4 exactly, because the product Plo(1−d/2)(1−εJ)(1−3d/4) is strictly larger than 1/2 even at the largest registered d. The code's integer expression was checked against an independent direct Fraction construction at 48 valid small test points; `verify_top64.py` additionally reconstructs every one of the 64 retained full parameter certificates and every worker's 32 certificates directly from the original scalar formulas. No displayed decimal is used to reconstruct a proof value.

Each 512-k chunk retains 32 candidates before the exact global merge. This normally would be insufficient by itself to claim a complete global top 64. Here a data-dependent check closes that specific gap: no chunk contributes 32 or more entries to the final list (the maximum is 21). If a discarded candidate ranked within the global top 64, its own chunk's 32 higher-ranked retained entries would all also rank in the global top 64, contradicting the observed occupancy. Therefore, conditional on the recorded frozen exhaustive enumeration and its exact comparator, the 64 listed values are the actual top 64 for this finite list. The verifier does not independently rescan all 7.37 million tuples or establish the conditional operator theorem.

The run used 15 distinct worker PIDs, 0.443 s elapsed and about 4.199 s total recorded worker CPU. The small registered warm slice was completed and reported before the full count. The entire exact enumeration finished before the first 0.5-second monitor sample, so `monitor.jsonl` is empty; no actual instantaneous CPU or memory sample was obtained and no full-CPU assertion is made. Worker-level cumulative CPU, RSS, per-worker top32 and chunk counts are in `results.json`. Registration, code, Fraction-equivalence output and resulting hashes are separately recorded.

Check the retained candidates from the clone root without repeating the full enumeration:

```
PYTHONDONTWRITEBYTECODE=1 /Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python cube-dim-7-bounds/upper/round052/compute/verify_top64.py
```

The verifier is a retained-candidate check and does not repeat the full search. The 7.37-million-tuple run itself is not repeated by this command.
