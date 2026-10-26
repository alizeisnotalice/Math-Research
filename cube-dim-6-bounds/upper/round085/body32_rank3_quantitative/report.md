# A085 report: a real body32 rank3 quantitative package at a new point

The package completes seven fixed local configurations at
x*=(1,1,1,1,1/2,0), with no old061 point, native routine,6D full kernel,
or32-vertex baseline. It supplies the actual selected-U amplitude through
C4, a uniform C2 time-remainder package, and exact quotient/basis data for
accepted076. **No usable complete Green semiconcavity constant Gamma was
obtained.** The negative-profile constant alone is not a Gamma for G.

The complete derivation and interface are in `quantitative_package.md`;
all exact tables are in `results.json`. The parent's final proof review is
still authoritative. The primary unique-U appendix was read and pinned
without rerunning its verifier; it awaits the independent review arranged
by the parent. Original registration and all three frozen source/input
files remain unchanged.

## Main exact results

The selected incident span is U={(a,a,a,a,b,c)}. In orthonormal tangential
coordinates y=(sumfirst4/2,x5,x6), its four original columns are(2,+/-1,+/-1),
and y*=(2,1/2,0). The complete true3D density satisfies

    pi³ p_t^U(y*)/t³=2560/[(1+16t²)(1+64t²)(9+64t²)],
    a_U(y*)=2560/(9pi³),
    pi³ r_U(0,y*)=-2007040/81,
    p_t^U=t³a_U+t⁵r_U.

The one-dimensional four-Cauchy fiber has Jacobian1/8; all four upper-half-
plane residues were included. The general amplitude formula

    a_U(y)=1/(8pi³) sum_j product_(ell!=j)(zeta_j(y)-zeta_ell(y))^-2

is evaluated with its full four-column normalization. Its35 partials through
order4 and35 uniform box bounds are stored exactly. Ten additional uniform
partial bounds through order2 control r_U on the full registered time/space
range, using a positive power-series majorant, not a sample grid.

The registered scales are

| Parameter | Exact value |
|---|---:|
| Tangential box center | (2,1/2,0) |
| Tangential box halfwidth | 1/32 |
| Time upper bound T | 1/64 |
| Local ball radius rho | 1/512 |
| Fixed scalar cutoff R0 | 1/41472 |
| Minimum pole-center separation on the box | 7/32 |
| True S2 distance squared at x* | 1/8 |
| Lower squared distance to a nonincident rank3 span | 1/864 |

The S2 result is an exact32-partition scalar two-level relaxation with an
actual body-pair witness, not an all-pair search. The nonincident distance
uses integer Gram determinants and Hadamard's bound, without enumerating
triples. The parent's independently sealed appendix proves uniqueness of U
at this point; that dependency is explicit in `external_dependencies.json`.

In the fixed orthonormal normal coordinates, the actual outside28 columns
project to three axis directions2e_i and four tetrahedral directions, each
with multiplicity4. Their marginal exponent is exactly

    psi(xi)=8||xi||1+max(16||xi||infinity,8||xi||1),
    16sqrt2 |xi|<=psi(xi)<=16sqrt3 |xi|.

The constants are sharp. Hence the deficit norm satisfies

    |z|/(786432pi)<=N_psi(z)<=3|z|/(1048576pi).

At y*, the actual selected negative profile a_U(y*)N_psi(z) lies between
5|z|/(13824pi⁴) and5|z|/(6144pi⁴). Only this marginal norm law is grouped.
All28 tangent/weighted labels remain separate, with an explicit basis
containing each label. Their summed coarea bounds are

    M0<=t^-3/(8pi³),  M1<=11t^-2/(4pi³),  M2<=11t^-1/(2pi³).

Seven distinct inverse matrices and56 exact dual-sign values certify the
norms;28 label mappings identify actual body columns. The fixed derivative
basis has quotient matrix2I3 and tangential transport operator norm1.

## Budget interpretation and remaining work

For f_U=-a_U(y)N_psi(z) alone, on the local ball, the rigorous upper
semiconcavity constant is

    3465893065395297320960 / [65143564997786663569923 pi⁴].

Its orientation-only decimal is0.0005461905056499521; the corresponding
Lipschitz bound is approximately1.4964826464204757e-5. Exact rational-times-
pi formulas, and a rational upper bound using pi>=3, are stored in Q007.

This result does not include the positive curvature of the actual Green
remainder. The full h4 extension C4 and correction C2 bounds are also supplied,
but are much larger than their natural-box counterparts. They are valid
existence-level data, not evidence that the assembled bound is competitive.
No parameter expansion was launched to conceal this limitation.

The unpaid ledger is the selected-U cutoff-removal, weighted quadratic and
finite-time signed-moment totals; the independently owned086 rank>=4
aggregate contribution; the original full-G t>=T tail; and compatible true
scalar vertices. The signed layer must retain its linear renormalization
and upper tail. The086 aggregate proof was not duplicated and no2^32-pattern
enumeration occurred. Conditional on the primary unique-U appendix passing
review, no additional incident amplitude needs to be explored for this point.

## Actual execution and provenance

One scientific program ran once, with seven configurations started and seven
completed. Guard PID91773, child PID91774, tool chunk e99b56. Checker CPU
0.233814 seconds and wall0.338497917 seconds; whole guard wall0.402595291
seconds. OS child user+system CPU was0.293191 seconds. Including guard CPU,
recorded program CPU was0.348957 seconds, below the30-second budget.

The guard had a180-second whole-program hard alarm,165-second child wait,
and10-second own-process kill/wait allowance. Its soft2/hard29 CPU setting
allowed the child to lower its inherited limit to27 seconds. No timeout,
startup failure, scientific failure, live child, or automatic retry occurred.

Actual mathematical work was:

* 32 column transforms and32 new-point scalar two-cluster assignments;
* one complete3D density path with four residues;
* 35 exact amplitude point partials,35 uniform amplitude partial bounds,
  and10 uniform time-remainder partial bounds;
* seven quotient inversions,56 dual-sign norm values, and28 label mappings;
* the registered cutoff, coercivity, and single-profile bound arithmetic.

There were zero6D/body32 full kernels, native calls, old baselines, MATLAB
calls, recursive agents, automatic retries, staging operations, commits,
or pushes. The frozen input, checker and guard hashes were checked before
execution and again at sealing. The primary appendix arrived after the
seven-configuration registration and was bound read-only in a separate
dependency record, as anticipated in the original input. No original pin
or scientific input was silently changed.

No whole-box certificate, continuous-I integral, global g bound, or Cn
improvement is claimed. The fixed data package is ready for independent
review and for a separately budgeted complete-remainder assembly.
