# r030 fixed-H receiver validation

This is a floating, finite-domain validation of the archived r029 `theta=0`
input; it performs no SQP or other reoptimization.  The code rebuilds the
receiver kernels from a replaceable JSON seed and splits every inner t2
integration interval at roots of the likelihood equations
`S=1`, `S=n`, and `S=Phi'(A)`.  These are respectively the H jump, the
start of the soft branch, and the height cap boundary.

With the original finite domain and 2049-point root scan, the root-split
results are:

| outer nodes | inner order | scaled cost | actual log J | minimum average |
|---:|---:|---:|---:|---:|
| 128 | 16 | 0.06078651151765332 | -68.48331463791287 | 0.998090661028606 |
| 256 | 16 | 0.06078651151670295 | -68.48331463792850 | 0.998090661026112 |
| 512 | 16 | 0.06078651151543926 | -68.48331463794929 | 0.998090661025706 |
| 256 | 32 | 0.06078651577413020 | -68.48331456788949 | 0.998090661033048 |
| 256 | 64 | 0.06078651577413017 | -68.48331456788949 | 0.998090661033048 |

`scaled cost` is `exp(tau_ref)*J`, while `actual log J` is
`log(scaled cost)-tau_ref`.  The 32-to-64 inner-order change is below
`3e-17` in scaled cost and `2e-16` in the minimum average.  Increasing the
inner order from 16 to 32 changes the minimum average by only `6.94e-12`.

The old unsplit 256-by-256 quadrature is separately recomputed by the
script and exactly reproduces the archived minimum 1.0006102883048078.
The root-split calculation instead gives about 0.99809066103, so the
archived strict feasibility margin is not robust to resolving the
H=0-to-1 discontinuity.  This is numerical evidence only; no tail bound,
interval arithmetic, or continuum feasibility certificate is claimed.
