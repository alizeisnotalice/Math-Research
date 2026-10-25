# Round 008: two-group tail-spacing break

The registered family fixes `a1=1+u^3/2` and uses

```
a2=a3=s+delta,   a4=a5=s-delta,   s=(1+u)/2.
```

Only positive `delta` was needed: changing its sign merely permutes the two
tail-coordinate pairs.  For every candidate, all 1024 length vectors in
`{1,2,3,4}^5` were processed in `Q(u,v)`, `u^5=2`, `v^5=3`.  Recursive
slicing and an independent coordinate-mask union agreed.  The comparison was
made against C's complete 25-coefficient algebraic value, not its displayed
decimal or an isolating-interval endpoint.

| delta | display only | exact conclusion | retained / nodes / seconds |
| --- | ---: | --- | ---: |
| 1/512 | 1.9909005901719705 | below C (upper difference < -0.00341587139651525) | 42 / 551 / 0.667 |
| 1/256 | 1.9876136043452675 | below C (upper difference < -0.006702857223218385) | 42 / 551 / 0.657 |
| 1/128 | 1.9814274051476077 | below C (upper difference < -0.01288905642087808) | 42 / 542 / 0.656 |

All jobs stayed well below the 120-second and 200,000-node checkpoints.  This
finite batch therefore gives no lower-bound improvement and needs no new
general-input or all-scale transfer interface.  It does not rerun the prior
common-tail, length-five, or alternating-weight baselines.

The next distinct mechanism should first find a non-dominated new witness
under a genuinely two-parameter tail split, such as `a2=a3=s+p` and
`a4=a5=s+q` with independently chosen rational `(p,q)`, before attempting a
two-parameter event-cell certification.  A symmetric opposite split alone is
now diagnostically unfavorable at the registered scales.

## Hashes

* registration: `4fed7e88db52212a7462671599c83a83b82bab1d69d15fae1f8f514285245a77`
* source: `32aad097557126865bdae072483e645c401371c5f7ca5caf75559125e1835dd8`
* delta 1/512 result: `a6b97164a03e1142ebbba33136a8a76ecbb36e52bea54d2ab40f1215d3847f54`
* delta 1/256 result: `2a6fe4747639a006a36f32cea9d35f3fdd36318baa3011506d7652e5d649b5a2`
* delta 1/128 result: `8fff261774f1e7b9f3762bb0d7683df2aff158255af4728840025106218f81f5`
