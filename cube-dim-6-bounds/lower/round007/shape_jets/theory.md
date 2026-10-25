# Whole finite-union first variation

This is an original finite-geometry argument, not an imported operator
theorem. The operator, source normalization and finite-L1 transfer remain
those of the six-dimensional project.

## Separated endpoint jets

Consider finitely many axis-aligned boxes B_k(t), in a fixed bounded
rectangle for 0<=t<t0, whose lower/upper endpoints have right expansions
`x(t)=x0+t*x1+o(t)`. First assume all pairs of distinct endpoints in a
coordinate have different pairs `(x0,x1)`, except endpoints already known
to be identical functions. Then lexicographic order of `(x0,x1)` is the
true endpoint order for all sufficiently small t>0: distinct constants
dominate their O(t) perturbations; equal constants with distinct first
coefficients are ordered by the first coefficient. Finiteness permits
one common positive t0. Thus all nonempty slabs and their active box sets
are fixed on this punctured interval.

Slice the union by the complete endpoint set in the first coordinate.
On each open slab the active box family is constant, so its volume is
slab width times the volume of the projected active union. Recurse in
the remaining coordinates and use ordinary interval union in dimension
one. This is an exact finite disjoint slab partition; all boundary faces
are null. Replacing addition/multiplication by

```
(a0,a1)+(b0,b1)=(a0+b0,a1+b1),
(a0,a1)*(b0,b1)=(a0*b0,a0*b1+a1*b0)
```

computes precisely the constant and right-linear volume coefficients.
Zero baseline widths with positive first coefficient must be included.
Two simultaneous thin directions give no first-order volume, by the
displayed product. Baseline contained boxes may cease being contained,
so the computation includes the full registered raw family.

## Higher-order ties: volume derivative versus true endpoint order

For general C1-right endpoints define their affine approximations
`xhat(t)=x0+t*x1`. Distinct true endpoints can have identical jets and
different second or higher coefficients. Lexicographic ordering does
**not** resolve their true order, and no equality for all t is asserted.

Nevertheless the volume coefficient is determined by these jets alone.
To see this, all original and approximating boxes lie in [-R,R]^n for
small t. For intervals, including empty ones,

```
|[l,u] symmetric_difference [lhat,uhat]|
 <= |l-lhat|+|u-uhat|.
```

The symmetric difference of two products is contained in the union of
the coordinate interval symmetric differences times the other bounded
coordinate ranges. Thus the volume difference of a box and its affine
approximation is at most
`(2R)^(n-1)*sum_coordinates(|l-lhat|+|u-uhat|)=o(t)`.
The symmetric difference of two finite unions is contained in the union
of corresponding box symmetric differences. Therefore the complete
union volume and the affine-box union volume differ by o(t).

For the affine boxes, lexicographic sorting is exact for small t, with
identical jets now being identical **affine approximations**. Their union
volume is a polynomial on the resulting eventual combinatorial region;
the dual-number recursion gives its first two coefficients. By the o(t)
bound these are also the true constant and right derivative. An interval
whose two affine endpoints coincide can be discarded from the affine
union, even if its true width is a positive o(t); its volume contribution
is o(t). This does not discard any first-order newborn slab.

Consequently higher-order ties leave true endpoint order and quantitative
finite step size unresolved, but do not leave the first derivative
unresolved. This extension is being separately checked by E. A positive
derivative still requires a new finite-step exact inward-volume certificate
before changing the C6 lower bound.

## Application to the registered motif

All selected masses are positive at t=0. Their sixth-root radii are smooth
near zero. All neighboring source gaps are positive at t=0, so consecutive
index endpoints remain the spatial extrema for small t. The source-endpoint
functions are polynomials of degree at most two. Clipping by the moving
cell faces uses finite min/max; their right expansions are obtained by
lexicographic min/max. Hence the hypotheses above hold uniformly over the
finite family.

For active start q, length l and initial radius r0,

```
r1=r0*D*(-1)^q/(6*l)    if l is odd, otherwise 0,
p_q,0=a_i*q,
p_q,1=a_i*(q*h_i+E*(q mod 2)).
```

Inactive endpoint derivatives are a_j times index times h_j. The active
cell upper-face derivative is 2*a_i*h_i; each inactive upper face has
derivative a_j*h_j/2. These are all in the old exact field. The ratio
has mass denominator 2 and reflection factor 32 at every t, so its
right derivative is exactly 16 times the computed union derivative.
