# r032 receiver-dependent continuous check

Status: floating-point numerical evidence only; it is neither an interval
certificate nor a claim about the unrestricted endpoint constant.

The source profile is the r031 repaired profile, with common source parameter
`c0=1.13434837109461`, fixed height cap 2, and split roots of the source
likelihood at `1`, `n`, and `Phi'(2)`.  Each root is obtained by monotone
bisection in `t2`.  The `H=1` and `H=2` pieces use Gamma--binomial CDF
differences and exponential cost primitives; only the soft piece uses
Gauss--Legendre integration.

The receiver-dependent test raises only receiver groups 23, 24, 25 from
`c0` to `c0+0.05`; the source likelihood, profile geometry, and all source
multipliers are initially untouched.  This is a direct numerical perturbation,
not a monotonicity argument: in `Gamma_b+c Binomial_b`, the Bernoulli
probability is `(1+exp(c))^-1`, which decreases with `c`, while the nonzero
Bernoulli shift increases.  The resulting distributions cross, so no
parameterwise first-order stochastic-order or coupling claim is valid.

At the r031 safe source scale, this moves the bottleneck from group 24 to
group 26 and gives N=128, soft-32 minimum average `1.0000832579286225` with
unchanged scaled cost `0.06361763303012608`.  The source scale was then
lowered once to `1.0471293307430003`.  The final split integrations were:

| outer / soft | scaled cost | actual log J | minimum average |
| --- | ---: | ---: | ---: |
| 256 / 32 | 0.06352833805958752 | -68.4391964772287 | 1.0000234626275157 |
| 256 / 64 | 0.06352833805958752 | -68.4391964772287 | 1.0000234626275153 |

Against r031's corresponding `0.06361763302943355`, this is a `0.14036%`
scaled-cost reduction, or an actual-log-J improvement of about `0.00140461`.
The calculation deliberately does not re-optimize the 33 source multipliers.
The two fixed-c mass-transfer probes in the JSON increase cost, so they are
only a diagnostic interface, not the result claimed here.

## Necessary global-c control

This local change is **not** evidence that receiver dependence is necessary.
With the same source profile and repaired r031 scale, a direct 256/soft-32
control replaces *every* receiver parameter by `c0+0.05`.  Its cost is
identical before source rescaling (`0.06361763302943355`), but its minimum
average is `1.0027117762328754`, versus `1.0000832579272974` for changing
only groups 23--25.  Thus the uniform receiver shift dominates this local
step on the tested scale.  The exact control output is
`receiver_continuous32_allc_compare.json`; no additional scale search was
performed.  For every nonzero group `j`, the group parameter applies to both
mirror receiver points `u[j]` and `u[j+32]`; the local test changed both
points for each of groups 23, 24, 25.
