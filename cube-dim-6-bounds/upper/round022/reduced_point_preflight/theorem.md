# Structural preflight at x0=(1,1,1,9/10,27/50,0)

## The point is off S2 and, more strongly, off S3

For the body32 sign-vector frame, S_r denotes the union of spans of at
most r columns. S2 contains only vectors with at most two coordinate
absolute values. The four values 1,9/10,27/50,0 exclude S2.

There is also a short exact S3 obstruction. If y is a combination of
three sign vectors, the signs of a given coordinate form one of eight
rows in {+/-1}^3. Modulo global negation there are four row classes,
represented by (1,+/-1,+/-1). These four representatives satisfy a
linear relation with coefficients all +/-1. Thus, for any four selected
coordinates, either two have equal absolute value or their signed sum
with some four +/-1 coefficients is zero. The assertion also applies
to spans of fewer than three columns by adjoining zero coefficients.

Select coordinates 1,4,5,6 of x0. Their values are (50,45,27,0)/50.
They are pairwise distinct in absolute value. Up to overall sign, their
possible signed sums have nonzero absolute values 22/50,32/50,68/50,
122/50. Hence x0 is off S3.

This argument gives a useful exact separation bound without enumerating
all three-column subspaces. S3 is contained in the union of the selected
coordinate hyperplanes y_i+/-y_j=0 and the selected signed-sum hyperplanes.
At x0 the minimum absolute pair-form value is 1/10, with normal norm
sqrt(2); the minimum absolute four-term signed sum is 11/25, with normal
norm 2. Therefore

    dist(x0,S3)^2 >= min((1/10)^2/2,(11/25)^2/4)=1/200.

This is a lower bound, not a claimed exact nearest-subspace distance.
In particular the known locally uniform four-large-parameter argument
applies: p_t(x)=O(t^2) locally near x0 as t decreases to zero, and K_T
is C2 there. These are regularity statements, not new kernel values.

## Exact group reduction

The signed-permutation group has order 2^6*6!=46080. A transformation
fixing x0 may permute the three equal positive coordinates and may
change the sign of the zero coordinate. Every other coordinate and
nonzero sign is fixed. Thus the stabilizer has order 3!*2=12 and the
point orbit has 3840 elements. The actual enumeration verifies this
and retains the multiplicity 12 for every distinct image.

Consequently 3840*21=80640 orbit/cone pairs account for all 967680
logical group/cone terms. Each has alpha=|det R|/product q_j and exact
nodes b_j=x.r_j/q_j. Sorting nodes loses no information because the
cone density and divided differences are symmetric in those nodes.
For the real density, globally negating all six nodes leaves the term
unchanged. We may therefore also merge each sorted tuple with its
negative, provided both full rational weights are retained. Self-paired
tuples must be counted only once. The checker verifies partner weights
and counts before this merge.

The compressed tuples are an exact arithmetic representation of the
full cone sum, not evaluations of p, G, H, or K. In particular, zero
and repeated nodes are retained and counted. A generic-only evaluator
cannot be silently applied when the zero-node count is nonzero.

## Candidate density denominator and the cancellation gap

With z=T^2 and b=(b1,...,b6), the real density term is

    alpha [z^3-e2(b)z^2+e4(b)z-e6(b)] /
                        product_(j=1)^6(z+b_j^2).

The full p_T is the sum of these terms divided by (2pi)^6. For each
exact nonnegative beta=b_j^2 let m_beta be its largest multiplicity in
any tuple. Then

    D_raw(z)=product_beta(z+beta)^m_beta

is a valid common denominator. Its degree sum_beta m_beta is a raw
upper bound. The prefactor (2pi)^-6 does not affect the denominator.
Zero beta contributes a factor z and must not be deleted merely
because the physical point is nonsingular. Repeated beta may arise
from opposite b values even when the original b nodes are distinct.

To obtain a reduced rational function, a later authorized stage must:

1. form each even numerator exactly from its signed node tuple;
2. retain all repeated and zero factors while lifting to a common
   denominator, or equivalently compute exact partial fractions;
3. sum the full rational weights with exact arithmetic;
4. compute and divide the exact numerator/denominator gcd, or prove
   the vanishing of every pole coefficient claimed to cancel;
5. verify the reduced identity against the original exact rational
   circuit before making root-count or derivative-sign claims.

Even an obvious within-term cancellation is a separate exact algebraic
step; the present raw bound does not assume it. Finite G regularity
implies the appropriate full principal parts must cancel, but gives
neither a computed reduced denominator nor its degree.

The structural audit stops before all of these cancellation steps.
If D_raw exceeds the earlier 512-degree budget, the reduced-degree
gate remains unresolved: exact cancellation would first need a separate
bounded release. No Sturm computation, new time, or Green pass follows
automatically. If the raw degree is within 512, that only bounds the
size of a future exact symbolic task; it is not an all-time certificate.
