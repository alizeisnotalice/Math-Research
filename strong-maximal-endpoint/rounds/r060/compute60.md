# r060 split-branch fixed-c reoptimization

Floating-point evidence only: no interval enclosure, no certified tail bound, and no claim of global optimization.

Fixed `N=64`, `b=32`, `cap=2`, `q=1.00001`, `c=1.236066745079915`, `tau=66.43897513452868`. For each gamma, 33 receivers use `v=0,...,16`, `d=gamma v^2/64`, and symmetric `(d+v,d-v)` points.

The evaluator splits every t2 fibre at `S=1`, `S=N`, and `S=Phi'(cap)`; constant branches use Gamma--binomial CDF differences and only the soft branch uses GL16. The `256` and `512` rows below preserve, rather than reoptimize, the selected 128-node weights.

| gamma | log J (128) | min average | score logV-logJ | gap to 71.79 | optimizer |
|---:|---:|---:|---:|---:|---|
| 0.2 | -68.5751979330 | 0.9999974616 | 71.1655428645 | 0.6244571355 | 0: CONVERGENCE: RELATIVE REDUCTION OF F <= FACTR*EPSMCH |
| 0.4 | -68.8008158867 | 0.9999991135 | 71.0996319941 | 0.6903680059 | 0: CONVERGENCE: RELATIVE REDUCTION OF F <= FACTR*EPSMCH |
| 0.6 | -68.8857333942 | 1.0000099065 | 70.9893604990 | 0.8006395010 | 0: CONVERGENCE: RELATIVE REDUCTION OF F <= FACTR*EPSMCH |
| 0.8 | -68.9134890557 | 1.0000098443 | 70.8750885263 | 0.9149114737 | 0: CONVERGENCE: RELATIVE REDUCTION OF F <= FACTR*EPSMCH |

Selected by the 128-node floating score: gamma `0.2`. Its exact finite-point cone-union formula (no spatial grid) gave volume `13.33437025346498`.

At 256 nodes the selected weights had minimum average `0.999997461593`. A common multiplier `1.00030756029938` was bisectioned against actual branch-split averages to reach `q`; the displayed 256 and 512 values use those corrected weights.
| preserved weights outer nodes | log J | min average | score | gap to 71.79 |
|---:|---:|---:|---:|---:|
| 256 | -68.574902899069 | 1.000010000000 | 71.165247830501 | 0.624752169499 |
| 512 | -68.574902899069 | 1.000010000000 | 71.165247830501 | 0.624752169499 |

A positive reported gap means the floating computation is strictly below 71.79 at the displayed quadrature; it is not a numerical certificate.

Final independent no-optimizer replay also uses 512 outer nodes and GL64 in the soft branch. The full gradient norm remains about 5.97e-4, and the projected KKT residual is about 2.76e-4 after scalar feasibility repair; convergence by relative objective decrease is not an optimality certificate.
