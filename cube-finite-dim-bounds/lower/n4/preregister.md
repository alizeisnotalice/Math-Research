# n4 batch 003 preregistration

Mechanism: retain the exact seed a=(1+sqrt(2)/2,(2+2^(3/4)3^(1/4))/4,(1+2^(1/4))/2,(1+2^(1/4))/2). Atom indexed k in Z^4 is at a*k+(0,(-1)^k1*q2,(-1)^k1*q3,(-1)^k1*q4)/2, with mass 1+(-1)^k1*w. Thus period cell has side lengths (2a1,a2,a3,a4), mass 2. This includes non-product position/mass coupling; it is not another spacing-only search.

Witnesses: rectangular index blocks with lengths 1..4 in every coordinate. Their bounding extrema and total mass determine a common radius r=mass^(1/4)/2 at threshold 1. Intersect receiver boxes with one full period cell; never inherit separate coordinate reflections. Compute full 4D union. Boundary equality is a limiting witness, requiring inward shrink before L1 transfer.

Finite resources / stop: first baseline closure and a signed q2,w grid q2/a2 in {0,.01,.03,.08,.15,.3,.5}, w in {0,.02,.08,.2,.4,.7}; then q=(q2,q2,q2) normalized by corresponding spacings with the same grid; at most 84 non-baseline cases. Terminate this batch after that fixed list, or after first >R4+1e-7 candidate for exact certification. If numerical union cost >30 seconds per evaluation, reduce to q/a in {0,.03,.15,.5}, w in {0,.08,.4,.7} and document. No claim of family exclusion from a failed finite batch. Python for union and exact seed, MATLAB for independent parameter/radius checks. Do not alter global STATE or ledger.

## Batch 003b (registered after 003 returned no candidate, before running 003b)

Address the rectangular-index-witness limitation using actual cube-captured atom subsets. Fixed 8 cases: single/diagonal normalized shift f in {.03,.15}, w in {0,.2}. For each case use scrambled Sobol 2^12 centers in the full cell (seed 2026102504); skip centers already covered by the rectangular union. Enumerate atom indices [-3,4]^4 and examine sorted distance prefixes up to 256 atoms; retain the strict-density-maximizing prefix if its average exceeds1. Add its receiver box to the existing union and recompute joint volume. This may discover nonrectangular index subsets. Stop after 8 cases, or a ratio >R4+1e-7 for certification. This finite sampling cannot exclude the full family.
