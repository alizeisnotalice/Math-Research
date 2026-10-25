# Round 006: one-axis alternating-mass finite periodic model

This round changes the input model, not just the dyadic precision.  On the
periodic lattice of period `2*a1 x a2 x ... x a5`, first-coordinate atoms
alternate between masses `1+theta` and `1-theta`; all other coordinate
directions retain unit weight.  Every period has mass two.  For a contiguous
block with first-axis start `j`, its exact mass is

`(l1 + theta*(-1)^j) * l2*l3*l4*l5` for odd `l1`, and `l1*l2*l3*l4*l5` for
even `l1`.  The program takes a strict inward dyadic fifth-root radius, so
`32*r^5 < block_mass` for every retained receiver box.

The receiver domain is `[0,2*a1] x product_{k=2..5}[0,ak/2]`.  The other four
coordinate reflections have factor 16 while the periodic input mass is two,
so the finite-model ratio is exactly `8 * union_volume` (not 32).  Starts
`j=-5,...,3` and all `l in {1,...,4}^5` were enumerated for each registered
theta.  Each run used exact `Fraction` decisions; recursive slicing and an
independent coordinate-mask union agreed.

| theta | raw / retained boxes | exact-ratio display only | exact comparison to C |
| --- | ---: | ---: | --- |
| 1/8 | 128 / 105 | 1.9302385343040613 | below, gap about 0.06407792726442436 |
| 1/4 | 116 / 106 | 1.8331202119206609 | below, gap about 0.161196249647825 |
| 1/2 | 125 / 115 | 1.6652894744605435 | below, gap about 0.3290269871079422 |

All three candidates completed in under 0.38 seconds and under 1,400
recursive nodes, versus their 120-second / 200,000-node caps.  The full exact
fractions, strict comparisons, mass/radius numerators, and provenance are in
the corresponding JSON files.  Hence this finite batch has no improvement.
It is not a weak `(1,1)` bound: a finite collar and smoothing/transfer proof
for this non-uniform periodic input would still be necessary even if it had
won.

The first launch encountered a Python negative-exponent type conversion before
any box was generated; it is recorded in `execution_log.md` and fixed by exact
integer parity before the three successful executions.

## Reproduction hashes

* `registration.json`: `9cf0c9cf4ee9b60715a5956cc64aca10671e4b51cc1457e78e2d1918c091e431`
* `alternating_first_axis.py`: `4154101294bd381dc98fc401f28f157bf27c1e6103de7d1919b8b827bac9bb15`
* `alternating_theta_1_8_results.json`: `2a900e39d1f4fe1f5df663b3f2ef8adc22d3e040a079e3c4543e4bbb6bfc3574`
* `alternating_theta_1_4_results.json`: `fa6bf4138789362fa9b786def838f4e3f847cb1c576aa558651ef50fe38ea26f`
* `alternating_theta_1_2_results.json`: `ec4e33d3242fa0c5adcd0ebaa59a8cb6e37d475d088386aee10b9a8dd2980f53`
