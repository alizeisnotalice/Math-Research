# Second-spacing candidate: finite L1 transfer

This applies to the candidate with

```
a=(1+u^4/2,(10001/10000)(2+u^5v)/4,(1+u)/2,(1+u)/2,(1+u)/2,(1+u)/2),
u^6=2, v^6=3, u,v>0.
```

For every retained block `ell`, put `m=prod ell` and
`r=m^(1/6)/2`.  The certificate's center box is the product of the exact
six slab intersections for the selected consecutive block.  Thus every
interior center has one common-radius cube containing at least `m` unit
atoms, while its cube volume is `(2r)^6=m`; possible extra lattice atoms
only raise its average.

Take compact inner boxes.  Finiteness gives a positive uniform boundary
margin.  Replace each atom by a nonnegative smooth unit-mass bump with
infinity-support smaller than that margin.  Every selected cube then
contains its full selected bump masses and belongs to `{M f>lambda}` for
each fixed `lambda<1`.

For an integer `K`, use receiver cells labelled `[-K,K]^6` and retain
source atoms labelled `[-K-2,K+2]^6`.  Every selected block has coordinate
indices in `[-2,2]`, so all its mass remains.  The source mass is
`(2K+5)^6`; the selected witness volume is
`(2K+1)^6(R'-o(1))`, where `R'` is the certified reflected-periodic union.
Hence

```
lambda |{M f_K>lambda}|/||f_K||_1
 >= lambda ((2K+1)/(2K+5))^6 (R'-o(1)).
```

Let `K` tend to infinity, the inner-box and bump radii tend to zero, and
`lambda` increase to one.  The conclusion preserves Lebesgue background,
centered axis-parallel cubes, one common radius, all scales, and finite
smooth `L1` inputs.  This is not a rectangle-scale statement.
