# Round 005: rational `l <= 5` finite-witness diagnostic

This is one new input, not a replay of the `l <= 4` baseline.  Its five
spacings are obtained by taking the exact lower isolating bounds of the locked
C spacings in `Q(u,v)`, with `u^5=2`, `v^5=3`, and rounding strictly downward
to denominator `2^40`.  Thus no displayed decimal is used to construct the
input.  For every length vector in `{1,...,5}^5`, the common radius is an
integer-computed dyadic lower bound `r=m/2^40` satisfying
`32*r^5 < prod(lengths)`; this puts every retained receiver box at strict
level above one.

The executed certificate considered all 3125 length vectors.  It found 56
nonempty boxes, five containing a length 5 coordinate, and after exact
dominance removal retained 42 boxes (none containing length 5).  Exact
`Fraction` recursive slicing and independent mask-grid cell enumeration agree
in 565 recursive calls.  It finished in 0.099 seconds, under the registered
120-second / 200,000-node limits.

The resulting exact ratio and the locked C upper endpoint are in
`rational_l5_results.json`.  Its display value is `1.9943164615594748`, which
is **below** C's exact upper endpoint by about `9.010926923587684e-12`.
Therefore this batch supplies no lower-bound improvement and does not make a
claim about the full weak `(1,1)` constant or about transfer.  The later
collar transfer audit remains separate; `±3` was only the registered safe
placeholder and was not invoked here.

## Reproduction and hashes

Run `python3 rational_l5.py` from this directory.  At the recorded execution:

* `registration.json`: `d660493ecda4c030df734437a410665e0189f3ae49bd9727e91dabc9dd0dc6ab`
* `rational_l5.py`: `e2af3878efd7ffd5e8af232b4057490660bac5e59b71fb8c89570c30bad8fb98`
* `rational_l5_results.json`: `47c85d2c3950dc48734e83e9bb84f56790af9bbc90c24fbec27a58bca7907b7e`

The runtime assertions verify the two exact union engines and every retained
strict radius inequality; the result includes the full rational ratio, C
endpoint, isolating-bound provenance, and all radius numerators.
