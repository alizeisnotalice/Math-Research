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


# Radius completeness plan for the fixed rationally rotated lattice

This is an append-only theoretical and computational checkpoint. No profile
job is launched until the parent hands off the serial lower-bound slot.

For cube side T, the three unrotated lattice counts are bounded by
(T/a+1)(T/s+1)^2. In the rotated plane, rescale by s and rotate back. A
captured integer point contributes a disjoint unit square centered there;
all such squares lie in the Minkowski sum of the rotated side-T/s square
with [-1/2,1/2]^2. The latter area is

    (T/s)^2 + 2*(7/5)*(T/s) + 1.

Thus the average of the full periodic atomic input is at most

    F(T)=(1+a/T)(1+s/T)^2
         *(1+2*(7/5)*s/T+s^2/T^2)/(a*s^4).

Every factor in the numerator decreases as T increases. An exact rational
verification of F(8)<1 therefore excludes every T>=8, uniformly in the
receiver. The script must actually verify this inequality before marking
the radius cutoff passed; this document does not substitute a numerical
display for that check.

For 4<=T<8, all changes in planar projection-endpoint ordering occur at
T=j*s/5, since every planar atom coordinate is an integer multiple of s/5.
The unrotated maximum count changes only at multiples of a or s. These
critical values, plus endpoints4 and8, give finitely many open intervals.
At an exact interior algebraic side, inner and outer rational arrangements
must have the same maximum planar overlap, with two independent algorithms
agreeing. That then fixes the planar maximum across the whole interval.
The maximum1D counts are also certified by separated rational intervals.

If the product of the three1D maximum counts and the planar maximum is at
most the fifth power of the lower interval endpoint, that entire open side
interval cannot witness a strict level1 average. Otherwise an upper cover
is formed as follows: enumerate permitted1D count triples, discard a triple
when its product times the planar maximum is too small, enlarge its canonical
receiver box to the upper side endpoint, and use the entire planar period
as an upper cover. Summing these product volumes is deliberately an upper
bound. It may suffice to show that unexamined larger cubes cannot close the
gap to C_cert without requiring a large new union calculation.

Critical side values themselves cause no omission for the strict level set:
any strictly winning cube at a critical side may be enlarged slightly while
preserving its strict inequality, and is then covered by the right-hand
open interval. The side8 and all larger sides are handled by F(8)<1.

For T<4, there are at most3 unrotated first-axis points and at most4 in
either unrotated tail axis. Those lengths are already in the registered
1..4 family. The only remaining small-side completeness question is whether
the planar count exceeds16. The actual exact maximum at side4 determines
the finite additional n values that must be checked. The 225 atom manifest
still suffices throughout: radius<=4 is smaller than the planar period5s,
so every atom in any relevant planar receiver square belongs to one of the
nine neighboring supercells.

These steps concern one fixed periodic input at one fixed angle. Even a
complete exclusion would not be an upper bound for the maximal operator on
general inputs or for other angles.


# Complete level1 enclosure for the fixed rotated input

The proposed fixed-input proof is now computationally complete and is supplied
for independent review. It concerns only the periodic unit-atom input

    P=(a Z) x (s R Z^2) x (s Z) x (s Z),
    a=1+u^3/2, s=(1+u)/2, u^5=2,
    R=[[4,-3],[3,4]]/5,

at threshold1. It does not bound the general maximal operator, other input
spacings, other angles, or other thresholds for this input.

For the centered axis-aligned cube maximal average over all positive radii,
the strict level1 superlevel volume per period, divided by period mass25,
has the rigorous rational enclosure whose decimal displays are

    [1.73995838416010979469265955500564837669736363224817715476779,
     1.73995838418300349717663336813054885976344568823264611332333].

The normalization is8/25 times the receiving volume in three unrotated half
periods and two full rotated supercell directions. No reflected rotated-axis
symmetry is used. The upper endpoint is strictly below the frozen complete
C_cert, so this one fixed input at level1 cannot be a new record.

## All larger cube sides are excluded

The exact Minkowski count argument in`radius_completeness.md` gives the
pointwise average bound F(T). Outward rational evaluation actually verified
F(8)<1. Monotonicity of every positive factor then excludes all T>=8.

The remaining range4<=T<8 was divided at every relevant j*s/5 planar
projection event and every relevant unrotated count event. The22 open
intervals exhaust this range up to their endpoints. In each interval,
complete endpoint ordering is fixed. An exact interior side received inner
and outer rational arrangements with the same maximum planar count;
independent bitmask-cell and event-sweep algorithms agree. The unrotated
maximum counts were separated by rational intervals as well. For all22
intervals, their product with the planar maximum is at most the fifth power
of the lower endpoint. Therefore no side in any open interval can have a
strict level1 average. A strict winning cube at a critical side could be
enlarged slightly without losing its strict inequality and would contradict
the right-hand interval. Thus all T>=4 are excluded, including critical
sides. This is a pointwise count argument, not only a sampled-radius test.

## Every smaller-side witness is already in the computed family

At side4 the planar maximum was actually certified to be15, by inner/outer
arrangements and both algorithms. Counts are monotone in side. Hence any
strict witness with T<4 captures at most15 planar points. Since a>7/4 and
s>1, its unrotated counts are at most(3,4,4), respectively.

For a one-dimensional lattice and receiver x in[0,b/2], the consecutive
length-l block closest in the maximum-distance sense is centered at0 for
odd l and at b/2 for even l. Its extreme-point distance is no greater than
that of any other length-l block. This is checked directly by comparing
the two neighboring possible block centers; all further centers are worse.
Thus replacing actual captured unrotated blocks by the canonical blocks
cannot destroy their common-radius inclusion.

The input is a product between those three lattice directions and the
rotated plane. If a cube captures counts(l1,l4,l5,n), its selected mass is
N=l1*l4*l5*n. A strict witness has T^5<N. Enlarge its side to N^(1/5): the
canonical unrotated blocks and at least n planar points are then included
at one common radius. The counts belong to the already executed lengths1..4
and n1..16 envelope. In fact N<=3*4*4*15=720<4^5, so the equality side also
remains below4. The225 neighboring atoms used by the program cover every
relevant planar capture region because the radius is below2<5s.

Consequently every strict level1 receiver belongs to the exact finite union
computed in v3. Conversely, each finite-family box interior can shrink its
equality radius slightly and has strict average above1. The finite union
of box faces has measure zero. The exact finite union volume therefore equals
the full strict superlevel volume for this input, up to null sets. Its
outward rational enclosure in v3 is thus a complete all-radius enclosure.

## Actual low-resource execution

Setup and22 interval checks ran strictly serially after the parent's explicit
handoff. Every job used nice10, one Python thread and six numerical-library
thread variables fixed to1. No NumPy or MATLAB process was used. Peak RSS
over these23 jobs was20,709,376 bytes, below20MiB. Each interval took roughly
0.007-0.010 seconds. Every result records its actual PID, nice value, thread
environment, peak RSS and wall time. The final saved-checkpoint verifier ran
as PID23012, nice10, one thread, peak RSS19,709,952 bytes and0.002231709
seconds. All these processes exited; there is no active lower computation.

The complete verification record is`complete_fixed_input_v4_results.json`.
Original K1/K8 and K2 checkpoints and their manifests remain unchanged. No
formal lower-bound promotion or Git commit was made. This closes the omitted
scale issue for this exact fixed input; it does not close the search for a
better lower-bound construction.
