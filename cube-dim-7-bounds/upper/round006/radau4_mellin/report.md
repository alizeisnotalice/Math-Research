# R006 four-moment Radau interval certificate

The resulting complete original-process integral budget is

    I < 7723627498344 / 10^18 = 0.000007723627498344.

The enclosed upper expression is approximately 7.72362749834377482e-6. It is strictly below a certified lower endpoint of the **R005 deterministic budget expression**, approximately 7.73042281751127158e-6; neither lower endpoint is claimed for the unknown actual I. This is about 0.088% below R005. It supplies no global g or C7 bound.

For the same tilted pair U and four-block target as R005, take W=cosh((log U1−log U2)/32)−1 and g(w)=T16(1+w)^(-1/4). Four exact moments of W use R004 q(1/8), R005 q(1/16), q(3/16), q(1/4), and exactly four newly certified q values:

| pair exponent r | q(r), approximate | new panels | Cauchy ratio h/ρ |
|---|---:|---:|---:|
| 1/32 | 0.891898204051347044 | 2048 | 2/3 |
| 3/32 | 0.710449776769789293 | 1024 | 4/9 |
| 5/32 | 0.566915508103060996 | 512 | 8/15 |
| 7/32 | 0.453159701568851226 | 512 | 8/21 |

The degree-two moment polynomial has two certified positive disjoint roots approximately 0.000975 and 0.005388, and all three corresponding Radau weights are strictly positive. Moment residual intervals through degree 4 contain zero. The frozen theory proves the true moment-defined Radau rule is exact through degree 4 and that the associated Hermite interpolant is a global upper majorant; the numerical ball polynomial alone is not that proof. Using lower endpoints of the two positive nodes and weights gives δ₆,L≈0.02430068458483039336 and therefore B₆≤1−δ₆,L. The original determinant and boundary interfaces yield M≤0.077148382596716926 and the I budget above. Full moment, root, weight, residual and dyadic correction endpoints are in certificate.json.

Each new q integral used positive [0,32] 512-bit acb_series, 129 coefficients, four initial pilot panels per exponent included in the 4096, an exponent-specific Cauchy remainder with the ratio shown, and the N=16 Gamma-recurrence bound for the complete tail. All 4096 unique panel records are terminal and positive. R004/R005 certificates are copied verbatim and their hashes checked against their source manifests; no old q panel was repeated. For the R005 comparison, the checker reconstructs an upper bound for its exact correction from upper μ/lower ν and combines it with lower q(1/8), giving a lower endpoint of that deterministic upper-budget expression.

Wall time was 27.920 s, with sum of panel CPU time 175.326 s, approximately 6.28 core-equivalents over wall. Five worker snapshots showed 15 processes summing to 639–662% CPU; measured summed worker RSS was 362256–498160 KiB, which may count shared pages repeatedly. The peak RSS of an individual process was 39,239,680 bytes; sampled system swap use was 1885.62 MiB. These observations do not claim system-wide peak memory or full CPU occupancy. Other work ran concurrently.

The checker validates hashes, 4096 unique positive terminal panels, independently reconstructs parameter-specific q error bounds and moments, verifies all root/weight and residual gates, then checks the correction and strict R005 comparison. It does not rerun the panel series coefficients; run.py is the reproducible generator. No additional panel, Monte Carlo run, MATLAB, commit or push was used.
