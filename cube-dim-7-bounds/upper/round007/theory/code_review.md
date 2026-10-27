# R007 read-only diagnostic implementation review

2026-10-27, n7_theory. Reviewed without execution:

    /Users/zhengzhihao/Desktop/T/Math-Research-cube7-matrix-r007/cube-dim-7-bounds/upper/round007/matrix_slack/run.py
    same-directory registration.md

The six proxy powers/factors, block identity detBg=2^18 product(w_i)sum(1/w_i), scaled full/pair Cholesky, and diamond EDGES match the theory. No total order between Zpairmat and Z1 is assumed. Global Y ratio flag and Cholesky diagonal flags trigger same-input128-decimal fallback. Original binary Y values are retained as64floathex strings and re-normalised at high precision; no clipping/redraw is present. A shared bounded semaphore acquired once per fallback and not released implements the global4096-row cap; exhaustion stops work.

The actual numerical registration explicitly chooses1024x1024-row chunks and SeedSequence([2026102707,stream_id]); these differ from the strategy's recommended chunk/seed example but were frozen by the numeric owner before the run and preserve the main-authorised2^20 rows perstream. Actual code retains6levels plus all15pairwise differences,21 observables, including signed unordered differences. Any simultaneous confidence reporting must use the actual list and the correct signed range2B where needed; the strategy's example list had15observables. This was communicated to the numerical owner.

All arithmetic/statistical output remains diagnostic, not a strict I expectation. This reviewer ran no checker, did not regenerate any shape and does not independently certify the completed sample evidence. Actual counts, fallback statuses, CPU and numerical values belong to Terra's receipts and the main audit.
