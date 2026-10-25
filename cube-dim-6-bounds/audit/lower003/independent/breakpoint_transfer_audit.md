# Breakpoint a2 candidate: audit and transfer

## Candidate and certificate status

The candidate changes only the second lattice spacing to

```
a2 = (u^3 + u^4 v^3)/4,    u^6=2, v^6=3,    u,v>0,
```

while retaining the other five spacings of the rectangular one-unit-atom
lattice family.  Hence it remains a Lebesgue-input lower-bound construction
for centered axis-parallel cubes with one common radius; it is not a change
to rectangle scales or the background measure.

The primary exact field certificate identifies it as an endpoint-order event
`r_(1,1,1,2,2,2)=2a2-r_(2,4,2,3,3,3)`.  Its all-candidate calculation has
127 nonempty boxes and 96 maximal boxes and encloses the raw reflected union
in

```
206692382601088 / 10^14 < R < 206692382601089 / 10^14.
```

That identifies the evaluated candidate; it does not establish global, or
even kink-local, optimality on the one-spacing line.

The separate main inward-rational checker uses rational lower/upper root
intervals, forms inward endpoints from the correct monotone interval ends,
and computes the full integer box union before optional containment removal.
For the `a2` expression every monomial is positive, so
`(u0^3+u0^4 v0^3)/4 <= a2 <= (u1^3+u1^4 v1^3)/4` is a valid enclosure.  Its
endpoint tests check both selected source-block ends against the same radius,
and `(2r)^6 <= m`; thus all witness cubes retain their common-radius
property.  Strict grid inward rounding provides a positive geometric margin.

The checker reads the frozen `rounds/003/k16_threshold.json`, not a mutable
ledger, and proves with exact fractions that both its inward union and its
literal finite-L1 ratio exceed the complete k16 exact upper endpoint.  Its
new strict finite target is

```
C6 > 5167309565027 / 2500000000000 = 2.0669238260108.
```

## Finite Lebesgue L1 passage

Put a unit-mass bump supported in the infinity ball of radius less than half
the recorded inward margin at every retained source lattice point.  For every
selected center, its same-radius witness cube contains every bump of the
selected consecutive block, and its volume is at most that block's total
unit mass.  Thus, for any fixed threshold below one, all selected inner-box
points are in the strict maximal superlevel set.

For receiver cell indices in `[-K,K]^6`, every selected block uses relative
coordinate indices in `[-2,2]^6`; retaining sources in
`[-K-2,K+2]^6` therefore contains the needed unit bumps.  The source mass is
`(2K+5)^6` and reflection gives selected receiver volume
`(2K+1)^6 R_in`.  The explicit finite ratio is

```
lambda * R_in * ((2K+1)/(2K+5))^6,
```

which the reviewed checker proves strictly exceeds the stated target (and
the complete k16 upper endpoint).  This supplies an ordinary finite `L1`
input directly, rather than relying only on an infinite atomic measure.

## Review boundary

The candidate lies at one enumerated endpoint-order equality.  Equality can
change the affine piece described in `k16_linearity_audit.md`, but it does
not invalidate the union or the inward certificate.  This audit does not
assert that the enumerated event list is exhaustive beyond the supplied
finite screen, or that the endpoint is the best possible spacing.
