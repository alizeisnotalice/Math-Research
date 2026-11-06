# r036 deformed two-block receiver curve

This is a finite-dimensional continuous split-quadrature model, not an
endpoint certificate and not a high-dimensional lift.  It changes the curve
geometry *and* rebuilds the source profile from reweighted multipliers.

Starting from the r033 scale (1.0276924765700282), common receiver value
(c_0+0.1075), and the r030 33-group weights, the two curve deformations are

\[
d_{\rm out}(v)=d(v)+0.15(v/16)^2,qquad
d_{\rm mid}(v)=d(v)-0.10v(16-v)/64.
\]

For each deformation, active groups 5--11 are multiplied by (e^\alpha),
while groups 19,20,21,33 are multiplied by
\(e^{-\alpha S_{\rm in}/S_{\rm out}}\).  Thus their combined multiplier
mass is balanced only to first order while the source likelihood and the KKT-profile geometry
are both reconstructed.  The exact Gamma--binomial CDF / exponential
constant-branch formulas and only soft Gauss--Legendre quadrature are used.

The outer-up bend is feasible and improves the tested cost; the middle-in
bend is infeasible in all three tested reweightings.  The selected tested
candidate is outer-up, \(\alpha=0.015\):

| quadrature | scaled cost | actual log J | min average |
| --- | ---: | ---: | ---: |
| 128 / soft-32 | 0.06223805808143694 | -68.4597158733417 | 1.0063184268973797 |
| 256 / soft-32 | 0.06223805807938099 | -68.45971587337472 | 1.0063184269092065 |
| 256 / soft-64 | 0.06223805807938098 | -68.45971587337472 | 1.0063184269092060 |

The undeformed counterpart has scaled cost 0.06239822910422077 and minimum
average 1.0029544314485879.

Each trial now retains only anchors whose numerical average exceeds
`1+1e-5`, then computes the analytic finite Pareto union of their positive
quadrants.  The baseline accepts all 33 groups, has volume
16.684838334873092 and score 71.27164607990836.  Outer-up also accepts all
groups but has volume 15.632170460599566 and score 71.20904687320880: it
loses 0.06259920669956 after volume settlement.  Middle-in accepts only
12, 11, 11 groups for the three tested reweightings, with respective scores
70.45270914960955, 70.41437500605355, 70.49266501013373.  Thus no tested
geometry wins the undeformed source.

By r026 equations (3)--(4), the deep-layer prefactor contains no
receiver-​(c) factor in these coordinates; no
\((1+e^{-c})^n\) adjustment is applicable.  The scalar reference is therefore
comparable, and outer-up remains 0.50594521474054 below its n=64 score
71.71499208794934.  These are numerical strictness checks, not interval
certificates or high-dimensional lifts.

On the fixed middle-in, \(\alpha=0.015\) source, three genuine
receiver-dependent directions were also evaluated with the mirror-swapped
coordinate rule (u_j:(c_1,c_2)), (u_{j+32}:(c_2,c_1)): the two signs of
\(c_1=c+0.1v/16,c_2=c-0.1v/16\), and the nonconstant common choice
\(c_1=c_2=c+0.1(v/16)^2\). Their accepted-cone scores are all
70.3760617997865, below the uniform-receiver value 70.49266501013373; their
minimum averages are lower as well. No refinement was run.
