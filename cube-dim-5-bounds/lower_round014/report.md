# Lower geometry 014: D5 full-parity finite screen

All six newly registered D5 candidates completed. Both union engines and the direct distinct-atom and ten-endpoint guards passed. Exact finite-family outer endpoints are strictly below the frozen complete C_cert; there is no improvement.

| Candidate | Exact inner ratio, decimal display | Nodes | Seconds |
|---|---:|---:|---:|
| anisotropic_27_32 | 1.69464812070471813650890500639757947163008107754327516455932 | 65192 | 1.444205500 |
| anisotropic_7_8 | 1.70179020148104991849473425341224109321533963766984175571922 | 29966 | 0.727125125 |
| anisotropic_29_32 | 1.58190216019035477818385410230666193052255925295329021804246 | 13100 | 0.409199750 |
| isotropic_31_32 | 1.56607028889922063652481957710889542984667232096515901947343 | 91289 | 1.320478459 |
| isotropic_1 | 1.59306072661983149054754161439095651131138250609119833630518 | 61404 | 0.981714000 |
| isotropic_33_32 | 1.56723337641206699410439656410265176598856091688723356374571 | 30131 | 0.516017042 |

Every exact equality-radius finite-family ratio is bounded above by its inner value plus `160*7776/2^80`. The exact rational upper comparison was executed. Full fractions and result hashes are in `completion_receipt.json`. The point set, mass formula, exact receiver intervals, reflection normalization and explicit width-one period collar are in `mechanism_and_scope.md`.

This excludes only the six spacing vectors and the chosen 7776 consecutive-block configurations. All other D5 parameters and omitted witness blocks remain unresolved. It is distinct from the sibling's two-axis parity at fixed C spacings, as confirmed before registration.

Reproduce from the repository root:

```sh
python3 cube-dim-5-bounds/lower_round014/d5_finite.py anisotropic_27_32
python3 cube-dim-5-bounds/lower_round014/d5_finite.py anisotropic_7_8
python3 cube-dim-5-bounds/lower_round014/d5_finite.py anisotropic_29_32
python3 cube-dim-5-bounds/lower_round014/d5_finite.py isotropic_31_32
python3 cube-dim-5-bounds/lower_round014/d5_finite.py isotropic_1
python3 cube-dim-5-bounds/lower_round014/d5_finite.py isotropic_33_32
```

Source SHA256: `c4fabbc543e3ad34ddd2d4d1e6dafa87fc1901a2d26a7c383e2736048d97ed25`. No commit or ledger update was made. No shared MATLAB process was used.
