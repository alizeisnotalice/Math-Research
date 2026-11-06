# r033 global receiver-​`c` calibration

Status: direct split-quadrature floating-point evidence only.  No
first-order-stochastic-dominance statement in `c` is assumed or used.

The source profile is held fixed at the repaired r031 source geometry and
source parameter `c0=1.13434837109461`.  Global receiver values were tested
at `c0+0.05`, `c0+0.10`, and `c0+0.20`, all at the r031 safe source scale.
The 128/soft-32 minimum averages were respectively:

| global offset | minimum average |
| ---: | ---: |
| 0.05 | 1.002711776235487 |
| 0.10 | 1.0037519701871176 |
| 0.20 | 1.0013684510412715 |

The only local refinement, offsets 0.075 and 0.125, gave 1.0034293251839017
and 1.0036940171642248.  Thus `c0+0.10` was selected.  One conservative
source-scale reduction (the only one in this round) from 1.0486657924183962
to 1.0276924765700282 was validated at outer-256/soft-32:

`scaled J = 0.06239822910366906`, `actual log J = -68.45714565701026`, and
`min average = 1.0029323661364373` (group 24).  This is a feasible-margin
probe, not a minimized scale or certified bound.

The originally archived asymmetric-block calculation was invalid: it used the
same parameter in both coordinate densities of a point, rather than assigning
one parameter per coordinate.  It is retained only under the explicit
`invalid_pre_correction_*` key in `receiver_global33.json` and must not be
used.

The corrected test is `receiver_global33_asymcorrected.json`.  For group 24,
the positive point is `u[24]=(11.904079694994362,-11.095920305005638)` and
its mirror is `u[56]=(-11.095920305005638,11.904079694994362)`.  On the first
point the density uses `(c1,c2)=(c_A,c_B)`, while the mirror uses
`(c1,c2)=(c_B,c_A)`.  This mirror-swapped mapping keeps the two receiver
blocks symmetric.  At the same reduced source scale and the same 128/soft-32
nodes, the `c_A=c_B=c0+0.10` baseline has minimum average
`1.002932366138988`; raising only group 24's `c_A` by 0.05 gives
`1.002739480861062`.  Thus the corrected one-round perturbation is negative;
no 256 rerun or additional search is warranted.  No source-multiplier
refinement or further geometry search was performed here.
