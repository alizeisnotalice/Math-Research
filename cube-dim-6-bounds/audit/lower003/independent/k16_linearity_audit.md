# k16 inward certificate and one-spacing breakpoint interface

## k16 certificate review

The reviewed parameter is

```
a_2 -> (10016/10000) a_2.
```

The reviewed main checker retains all relevant original features: rational
lower/upper isolating intervals for the sixth roots, strict inward endpoint
rounding at grid denominator `10^35`, endpoint inclusion checks against both
ends of every selected lattice block, and an integer recursive *union* rather
than a sum of overlapping box volumes.  It checks 127 inward nonempty boxes
and retains 96 after containment removal.  Its stated strict target is

```
258204069395 / 125000000000 = 2.06563255516.
```

The output makes exact `Fraction` comparisons against the full frozen R6
upper endpoint and the complete R5/R4 upper endpoints, and all are true.  A
dynamic-ledger issue is that this script alone does not compare against the
already certified k1 inward value.  `compare_k16_to_k1.py` corrects that
omission using the original exact rational results: both the k16 exact union
and its displayed finite-L1 target are strictly larger than the full k1
exact rational inward union.  It
also constructs a literal finite Lebesgue witness: uniform unit-mass bumps
of halfwidth `10^-35/2`, receiver indices `[-2^44,2^44]^6`, and source
indices `[-2^44-2,2^44+2]^6`.  The source collar contains every selected
relative block; the positive inward margin keeps the full bumps in the same
common-radius cubes.  Consequently its finite ratio is strictly above the
target.  This is a strict lower-bound audit, not an upper-bound result or an
optimality claim.

`replay_main_inward_k16.py` executes this checker read-only into the present
clone, redirecting only the result file and preserving the reviewed source
hash.  It does not alter the main clone.

## Why the one-spacing union is piecewise affine

Let `s=a_2` vary and hold the remaining five spacings fixed.  Do **not**
first remove contained boxes: use the union of every candidate box whose
five fixed-coordinate intervals are nonempty.  Once the endpoint branch in
coordinate 2 is fixed, each coordinate-2 endpoint has the form `alpha s +
beta`; all endpoints in coordinates other than 2 are constant.

Partition the five fixed coordinates by every constant box endpoint.  On
each relatively open resulting five-dimensional cell, the set of boxes that
can meet the cell is constant.  Its coordinate-2 section is a finite union
of intervals with affine endpoints.  If their endpoint order is fixed, the
merge pattern of this one-dimensional union is fixed, so its length is
affine in `s`.  Multiplying by the fixed five-dimensional cell volume and
summing over the finite partition proves that the full six-dimensional union
volume is affine on that interval.  Thus there is no product of two varying
lengths: only the coordinate-2 length varies.  This Fubini argument is the
needed correction to a generic ``six-dimensional box volume is polynomial''
objection.

The claim deliberately excludes breakpoints and says neither local nor
global optimality.  Using a changing retained-96 list alone is unsafe;
containment changes do not change the full union but can make an artificial
piecewise formula appear to jump.

## Exact breakpoint interface

For a count vector `ell` with mass `m=prod ell`, put `r=m^(1/6)/2`.  The
only coordinate-2 endpoint functions that can occur are

```
0,  s/2,  r-(ell_2-1)s/2  (ell_2 odd),
ell_2*s/2-r               (ell_2 even).
```

An exhaustive finite breakpoint enumerator must retain only vectors with
nonempty fixed-coordinate intervals and then include positive `s` values
from all of the following:

1. equality of the two alternatives in every `min` or `max` endpoint;
2. equality of the lower and upper endpoints of a candidate interval
   (zero-width birth/death);
3. equality of any two endpoint functions in the displayed finite set.

Items 1--2 are also covered if all pre-branch affine functions are included
in item 3, but listing them avoids accidentally omitting a clamp or a
zero-width event.  Sort these values by exact algebraic comparison, select
one rational/algebraic test value in each adjacent interval, and evaluate the
full union of all candidates there.  The affine formula then implies an
interval maximum occurs at a boundary (via one-sided limits); a strict lower
certificate still requires an interior parameter with positive endpoint and
smoothing margins.  This is a finite next-kink interface, not authorization
for unbounded stepping.
