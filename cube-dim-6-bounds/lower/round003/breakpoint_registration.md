# Round 003: second-spacing breakpoint batch

Only the second spacing of the `k=16` lattice varies.  The other five
spacings, unit masses, common-radius rule, and `ell in {1,2,3,4}^6` are
frozen.  For a fixed block, each second-coordinate endpoint is affine in
the spacing multiplier.  If all raw-box endpoint branches, zero-width
states, and pairwise endpoint orders are fixed, partitioning the other
five coordinates into their fixed cells makes the union a sum of affine
one-dimensional union lengths times fixed five-dimensional volumes.

The finite event list includes each min/max branch switch, every raw-box
second-coordinate zero-width event, and every equality of two possible
second-coordinate endpoints.  It does not use contained-box switches,
which do not change the full union.  At most 128 event values are allowed;
the first pass evaluates the nearest 16 around `k=16`, rebuilding the
whole union at each.  A result must clear the full `k=16` exact upper
endpoint by `10^-7` before one exact/inward certification is permitted.
