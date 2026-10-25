# Rational rotation of two equal tail axes: model fixed before execution

Let u=2^(1/5), a=1+u^3/2, s=(1+u)/2, exactly the frozen C_cert spacings.
Rotate only the equal second and third axes by

    R=(1/q) [[p,-r],[r,p]],  p^2+r^2=q^2,

using the single primitive triple (p,r,q)=(4,3,5). The other three axes have
spacings a,s,s and stay unrotated. Every atom has mass one. The orthogonal
matrix is invertible, so distinct five-integer index tuples remain distinct
atoms.

The rotated planar lattice is s R Z^2. Its axis-parallel supercell is
[0,q*s)^2, because (q*s,0)=s R(p,-r) and (0,q*s)=s R(r,p). The index of
these two integer vectors in Z^2 is det[[p,r],[-r,p]]=q^2. Equivalently,
write a planar lattice point as (s*U/q,s*V/q). Its original indices are

    m=(p*U+r*V)/q^2,  n=(-r*U+p*V)/q^2.

Thus representatives are precisely integer pairs U,V in [0,q^2) satisfying
both divisibility congruences. There are q^2=25 representatives, which the
source enumerates and verifies by exact integer arithmetic. They have
unique physical coordinates. Every planar lattice point can be uniquely
translated by multiples of q*s in each coordinate back into this half-open
cell; the inverse formulas then put it in the representative list. Hence
this is both a nonduplicated representation and a complete translation
cover. There are exactly 25 atoms per five-dimensional period cell, since
the other axes each contribute one atom.

Take the receiver domain, in computational order,

    D=[0,a/2] x [0,s/2] x [0,s/2] x [0,q*s] x [0,q*s].

Only the first three (unrotated) coordinates are reflected. Those eight
copies fill a five-dimensional supercell up to null faces, while the two
rotated coordinates are already used in full. If V is receiver union volume,
its finite periodic weak ratio is 8*V/25. Individual reflected rotated axes
are not assumed to preserve the input; no factor32 is used.

The finite family uses all lengths l_i in {1,...,4}. On each unrotated axis
choose the consecutive block centered about zero for odd length and about
a half-grid point for even length. On the rotated axes choose consecutive
original m,n blocks. Their starts are each of the 25 supercell representative
index pairs plus the nine supercell translations (h,k) in {-1,0,1}^2; the
corresponding start-index displacement is (p*h+r*k,-r*h+p*k). Every such
block has N=product(l_i) distinct unit atoms. Its one common equality radius
is r_N=N^(1/5)/2 in all coordinates. The receiver box is the product of
[max_selected coordinate-r_N,min_selected coordinate+r_N], clipped to D.

These nine translations cover every translated block of this registered
shape that can meet D: N<=1024 gives r_N<=2, and q*s>2. The start atom of
the representative lies in [0,q*s)^2 and belongs to its block. A translation
by at least two supercells in a coordinate puts that atom farther than two
from the whole receiving interval [0,q*s], so no common-radius witness box
can meet D. This checks translation coverage of the finite block shapes;
it is not an assertion that these shapes exhaust the full maximal set.

For actual computation, u,a,s and each N^(1/5) receive outward rational
brackets at denominator 2^80, verified using exact integer fifth powers.
Signed-coordinate extrema are propagated with outward interval endpoints.
The receiver boxes are rounded inward for lower witnesses and outward for
upper containment at a separate rational coordinate grid of denominator
2^40. The two unions therefore bracket the exact finite-family value for
the exact algebraic C spacings; no rounded spacing is substituted as input.
The inner radius is strict, so selected mass exceeds cube volume. Direct
integer-tuple enumeration checks distinct N and all ten receiver endpoint
inequalities for every retained lower witness. Two exact union engines are
required, and the entire candidate has a 120-second/200000-node budget.
