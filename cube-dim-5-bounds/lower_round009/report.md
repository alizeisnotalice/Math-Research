# Round 009: equal-mass alternating first-axis gaps

This is distinct from the unsuccessful alternating-*weight* experiment.  The
first-axis atoms all have mass one and, in one period, are at `0` and
`a1+delta`; equivalently `x_j=j*a1` for even `j` and `x_j=j*a1+delta` for odd
`j`.  The period length is `2*a1`, its total input mass is two, and the last
four axes retain spacing `s=(1+u)/2`.

All blocks `l in {1,...,4}^5` and first-axis starts `j=-6,...,5` were checked.
The receiver cell is `[0,2*a1] x product [0,s/2]`.  Tail reflection gives
factor 16 and the period mass is two, hence the exact finite-model
normalization is `8 * union_volume`; no cell-volume factor was introduced.
At equality radius the maximal average is at least one, and the standard
strict-level conclusion is obtained by taking thresholds below one and a
limit.  Both union engines are exact in `Q(u,v)` and agree.

| delta | display only | exact C comparison | raw / retained / nodes / seconds |
| --- | ---: | --- | ---: |
| 1/32 | 1.9907711433786968 | below (difference upper < -0.003545318189789112) | 129 / 109 / 821 / 1.962 |
| 1/8 | 1.9691275519898284 | below (difference upper < -0.025188909578657354) | 133 / 109 / 1038 / 2.246 |
| 1/4 | 1.8796291525229263 | below (difference upper < -0.11468730904555952) | 137 / 114 / 1146 / 2.608 |

All three finite candidates fail against the full 25-coefficient algebraic C
value, not merely its display or interval endpoint.  Thus this batch has no
new lower bound, no locking, and no transfer claim.  It stays below the
120-second / 200,000-node limits.

## Hashes

* registration: `2952773ea22ee09e682cd38d781a178cab8af80c1656232c5a249607c24650b1`
* source: `1efeef8cec90cab8372f53f4ad90dd5b2582d86325eb12a67f309ac434558cc6`
* delta 1/32: `92f0971d7b3f6ab92d3e0f91e70fc3f7c0da4a0bd4a200c7fd439739706ecd8a`
* delta 1/8: `e777489efc24b783f78985be49aca1dec5971d7790354d87e7a162aaace3a388`
* delta 1/4: `3774f9a676f98c81fc2e7751af185f2a320c9c2fe2de8c13e186b31edd113f1b`

## Next finite candidate

The tested gap modulation loses already at `delta=1/32`.  A different next
model should alter a short *two-axis* point pattern while preserving a finite
mass period, and first screen whether its new receiver boxes are non-dominated
before attempting exact unions.  Repeating this one-axis gap family on a
finer mesh would not be a distinct mechanism.
