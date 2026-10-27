# Proposed frozen pair-Mellin finite registration for Terra

Source theory: same-directory `pair_proof.md`. Scientific implementation/actual receipt belongs to Terra's `upper/round004/pair_mellin/`. This Sol strategy directory contains no numerical run.

Fixed model: original7D body64, equal positive weights, common Cauchy time; 8 disjoint tight Hadamard blocks. Target: certified q=E(T1+T2)^(-1/8), bound M<=q^4/2, complete I via105·2^(25/6)pi^(13/12) denominator. No g/C7 claim.

Parameters: s=1/8,a=1/16,T=32,N_tail=16,r_Cauchy=3/64,h_panel=1/32,512 equal panels[0,32],degree128 (series order>=129),512-bit outward complex arithmetic. Panelj centre=(2j+1)/32, j0..511. All constants exact rationals except explicitly enclosed Gamma/pi/log2.

Analytic majorant: M*=Gamma(1/64)Gamma(7/64)m(1/64)m(7/64).
Raw integrated Taylor remainder:32 M* (2/3)^129/(1-2/3).
q tail:m(1/16)^2 Gamma(257/16)^2/[pi Gamma(1/8)·31·32^31].
Assemble both errors as specified by pair_proof.md(9), then q→M→I(10).

Stages:
1. Interface probe: inspect available acb_series gamma/exp/inverse/order/precision APIs; no scientific panel claim.
2. Frozen4-panel pilot j0..3 with all coefficients and interval diagnostic identities. No MC and no old baseline node computation. Source/config/hash lock before scientific execution.
3. If pilot passes, fixed remaining508 panels distributed over15 processes, BLAS/OpenMP1, reuse pilot4. No wall-time cap; finite parameter/budget cap is the fixed panel set. Record true CPU/wall/peak concurrency. Do not claim15 busy workers without receipts.
4. Assemble all512 unique panels, verify exact coverage/indices/source hashes/degree/precision, rigorous full tail, exact outward intervals and strict comparison to R002 deterministic budget expressions. Round I upper endpoint to a specified rational grid only after certification; save its exact finer endpoint too.

Stop after registered assembly and main-agent review. If implementation fails or interval separation is inconclusive, report the exact issue; no automatic degree/tail/precision/sample expansion, no midpoint certificate, no submit/push by Terra. No Gamma numerical value or M_pair result is asserted in this unexecuted registration.
