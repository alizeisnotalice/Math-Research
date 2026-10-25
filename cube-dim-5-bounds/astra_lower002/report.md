# Astra lower002: full-displacement boundary kernel and one exact local obstruction

This round implements an analytic displacement response, rather than another
half-grid or small-mass grid. It does not improve C. It gives a strictly
negative first-order response on one explicitly bounded non-half-period
displacement interval, and leaves all other cells open.

For the exact standing lattice, add a satellite orbit of mass ε per period
at displacement z. The source period mass is1+ε. The normalized derivative is

    D(z)=B(z)+κ(z)−C,

where B integrates the largest active outward common-radius velocity over
the exposed union boundary, κ is the isolated satellite halo contribution,
and−C is the denominator derivative. The full-position formula and its
finite hyperplane arrangement are proved in `kernel_derivation.md`.

On every open displacement cell the response is multiaffine, total degree
at most four. This is proved from the ordering of the original receiver
endpoints and moving satellite-count cuts. It is not inferred from samples.
The actual exact run certifies the slice z=(0,t,0,0,0), restricted to the
complete registered arrangement interval containing a_2/4:

    L=(-v+u v²)/2 ≈0.26843575917524304,
    U=(1+u+u²−u³v³)/2 ≈0.2690251070932039,
    u⁵=2, v⁵=3.

There are1196 cuts in the safely over-refined slice arrangement, from145 raw
endpoint alternatives and36 radii. The three registered internal positions
are L+θ(U−L), θ=1/4,1/2,3/4. This is one finite cell calculation, not a search
over1196 positions. All these satellite positions lie strictly inside the
original unit-atom halo, so κ=0 on the whole interval.

The six reflected-halfcell kernel runs use exact Q(u,v) arithmetic. Each
takes about0.71--0.81 seconds and at most1692 weighted-union nodes. They
include all51 feasible or degenerate primary raw boxes;60 oriented face
groups survive trivial blockers. The derivative values are affine exactly,
including the full25-basis coefficient midpoint identity. The limiting
endpoint values are:

| Limit | Decimal display only | Exact conclusion |
|---|---:|---|
| D(L+) |−0.28741378179536053|Rational isolating upper endpoint<0|
| D(U−) |−0.2873531674634645|Rational isolating upper endpoint<0|

Therefore D(t)<0 throughout L<t<U for the defined finite witness family.
For each fixed t there is an unquantified positive ε range on which this
family remains below C. The same holds uniformly on compact subintervals.
The two endpoint strata, other positions and finite larger masses are not
included. This is not an upper bound on the full periodic maximal function.

An independently written mixed-receiver-box jet calculation checks the
middle node without importing the exposed-facet implementation. It rebuilds
14336 configurations over four nearest-block receiver intervals, yielding
396 raw jet boxes and162 retained boxes. Its constant term equals the full
C coefficient dictionary and its normalized derivative equals the facet
kernel coefficient dictionary exactly:

    D(midpoint) ≈−0.28738347462941255.

It verifies22678 distinct tagged atoms and6592 coordinate endpoint jet
inequalities, with noncoincidence of the two physical orbits proved by the
nonintegral displacement. Runtime was10.35 seconds and4430 union nodes.
Only the scalar algebraic field is shared with earlier code. The use of
lexicographic endpoint jets follows the idea in `lower_alt014`; this is a
separate same-agent implementation check, not an independent agent audit.

`kernel_derivation.md` also gives the finite collar, exact mass denominator,
inward receiver shrink and smooth-L¹ transfer for positive ε. All radii stay
common to the five coordinates. No record candidate is locked or promoted.

Reproduce from the independent clone root:

```sh
python3 -B cube-dim-5-bounds/astra_lower002/kernel.py register
# For i=0,1,2 and sign=1,-1:
python3 -B cube-dim-5-bounds/astra_lower002/kernel.py evaluate i sign
python3 -B cube-dim-5-bounds/astra_lower002/kernel.py finish
python3 -B cube-dim-5-bounds/astra_lower002/jet_check.py
python3 -B cube-dim-5-bounds/astra_lower002/verify_results.py
```

The next useful extension is a deliberately selected different displacement
cell, especially one where κ=1 or simultaneous multi-coordinate captures
change the active maximum. Neither this local obstruction nor the earlier
27 half-period diagnostics excludes such a mechanism. No additional cell
was sampled in this round. Only local Python was used; MATLAB and the main
clone were untouched. No commit or push was made.
