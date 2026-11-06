# r034 receiver-coordinate audit

All entries are direct Wolfram split-quadrature evaluations at the fixed r033
source scale `1.0276924765700282`, with outer 128 and soft order 32 unless
otherwise stated.  They are floating numerical evidence, not certificates;
no stochastic-order argument in `c` is used.

The correct group-24 mapping is `u[24]:(c1,c2)` and
`u[56]:(c2,c1)`.  At r033's common receiver value `c0+0.10`, the common
baseline minimum average is `1.002932366138988` (group 24).  Direct one-axis
coordinate probes give the following useful behavior:

| coordinate perturbation | minimum average | bottleneck |
| --- | ---: | ---: |
| c1 + 0.005 | 1.002940553706471 | group 25 |
| c1 + 0.010 | 1.002940553706471 | group 25 |
| c1 + 0.0125 | 1.002940553706471 | group 25 |
| c2 + 0.005 | 1.002938863372999 | group 24 |
| c1 - 0.005 | 1.002913273664705 | group 24 |
| c2 - 0.005 | 1.002917696655906 | group 24 |

Larger positive and negative steps through 0.05 decrease the minimum.  The
central direct-difference response at step 0.005 is positive in both
coordinates, but the minimum switches from group 24 to 25, so the min
objective is nonsmooth.  The best coordinate value was independently stable
at outer-256 with soft orders 32 and 64:
`1.0029405537044587` and `1.0029405537044584`.

However, this is not a receiver-specific gain.  Necessary uniform receiver-c
controls at offsets from r033's global value `c0+0.10` were tried at
`+0.005,+0.010,+0.0125,+0.025`; the best of those was +0.010 with
`1.0029542338453086`.  One additional bracket refinement +0.0075 improved
this to `1.0029544314485879`, still bottleneck group 24.  Thus the uniform
global receiver parameter beats the best local coordinate perturbation at
the same source and quadrature grid.  This rules out recording a distinct
receiver-coordinate mechanism from this fixed-source family; the next
discriminating direction is geometry rather than more fixed-source c fitting.
