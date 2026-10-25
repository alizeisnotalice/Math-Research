# Round 012: isolated two-atom clusters cannot approach the current bound

This is a new finite-cluster screen, not another spacing or weight grid on
the failed periodic families. The outcome is a structural exclusion for two
atoms of arbitrary positive masses, followed by three actual exact checks.
No new C5 lower bound is claimed. The proof below is a proposed elementary
argument for independent review by the main agent.

## Exact finite witness model

Let mu=w0 delta_p0+w1 delta_p1 with distinct p0,p1 in R^5 and positive
weights. At level alpha>0, each nonempty subset I of {0,1} has mass
m_I=sum_{i in I} wi, and one common centered cube radius

    r_I=(m_I/alpha)^(1/5)/2.

Its receiver box is the product of the five intervals

    [max_{i in I} p_iq-r_I, min_{i in I} p_iq+r_I].

Empty or zero-volume boxes contribute zero. There are just three subsets.
Up to null boundaries their union is exactly {M mu>alpha}, rather than
merely a selected lower witness family: any witnessing cube contains either
one atom or both; its side must be smaller than (captured mass/alpha)^(1/5).
Conversely, a receiver in the interior of a subset box can use a slightly
smaller common radius, giving strict average above alpha. Box boundaries
have zero Lebesgue measure. Thus this finite model handles all positive
scales for the two-atom input.

## General two-atom obstruction

The following argument works in every dimension n>=1, not only n=5. Scale
mass by alpha and write W=(w0+w1)/alpha, R=W^(1/n). Translate and reflect
coordinate axes so that the heavier atom is at 0 and the other at
(d1,...,dn), with di>=0. Let L=(max(w0,w1)/alpha)^(1/n), so L<=R and
L^n/W>=1/2. Let B be the two-atom receiver box, with side lengths
bi=(R-di)_+, and V=|B|<=R^n=W. The singleton receiver cubes have total
volume W, so their union has volume at most W.

If any di>=R/2, then V<=R^n/2=W/2 (including the empty-box case). The
whole three-box union consequently has volume at most 3W/2.

Otherwise every di<R/2. The intersection of B with the heavier atom's
singleton cube has coordinate length

    ci=min(L, (R+L)/2-di).

Both alternatives are at least (L/R)(R-di). The two exact residuals are

    L-(L/R)(R-di) = L*di/R >= 0,
    (R+L)/2-di-(L/R)(R-di)
        = (R-L)*(R-2di)/(2R) >= 0.

Multiplying positive coordinate inequalities gives

    |B intersect heavy_singleton| >= (L/R)^n V >= V/2.

The portion newly added by B to the singleton union is therefore at most
V/2<=W/2. Again the full union has volume at most 3W/2. In all cases,

    alpha * |{M mu>alpha}| / (w0+w1) <= 3/2.

This proves a ceiling for the entire two-atom atomic model, at every alpha,
for arbitrary displacement and positive masses. It is not an upper bound
for general L1 inputs, arbitrary smooth bumps, or configurations with three
or more atoms. The argument does not claim that 3/2 is sharp for two atoms.

## Consequence for separated groups and periodic models

For a finite input whose atoms are partitioned into one- or two-atom groups,
write Wtot for total mass. At a specified level alpha, assume that the
sup-norm distance between every pair of atoms belonging to different groups
is at least (Wtot/alpha)^(1/5). A cube with average strictly above alpha has
side less than (Wtot/alpha)^(1/5), so it cannot meet two groups. The full
level set is then the union of the groups' level sets. Subadditivity and the
two-atom theorem bound the total level ratio by 3/2.

For infinitely periodic groups there is no finite Wtot, so this finite
separation argument does not apply to arbitrary cross-group large cubes.
Still, any proposed periodic certificate using only within-group witnesses
has certified volume per cell at most the sum of the within-group volumes,
and hence ratio at most 3/2. Exceeding the standing C therefore requires
at least one of: genuinely interacting three-or-more-atom groups, or
cross-group witnesses at the target level. Enlarging the spacing of isolated
two-atom groups cannot by itself produce the desired improvement. This is
only a necessary condition; neither option is asserted to suffice.

## Registered finite checks and exact formulas

All examples have total mass one and alpha=1. The first atom is at 0.
The registration predates every run and contains no more than three inputs.
For every one of the three subsets, the program builds exact rational inner
and outer receiver boxes. Integer fifth-root bracketing gives side bounds
qlo^5<m<=qhi^5 at denominator 2^80. Thus inner boxes are strict lower
witnesses and outer boxes contain the full atomic level set. Direct selected
atom endpoint checks and two exact union engines (recursive sweep and
inclusion-exclusion) agree in each run. These finite checks corroborate the
implementation; the general exclusion follows from the argument above.

1. Equal masses 1/2,1/2, displacement (15/16,0,0,0,0). Put L=2^(-1/5).
   Since L<7/8, the two singleton cubes and the pair box are disjoint except
   possibly on boundaries. Their volumes sum to 1+1/16=17/16 exactly.
2. Equal masses 1/2,1/2, displacement (1/2,...,1/2). The pair volume is
   1/32. Its intersection with each singleton has volume L^5/32=1/64;
   the singleton-overlap volume equals the triple-overlap volume. Exact
   inclusion-exclusion therefore gives volume one.
3. Masses 31/32,1/32, displacement (15/16,0,0,0,0). Let L=(31/32)^(1/5).
   The singleton cubes are disjoint, the pair box has volume 1/16, and it
   meets only the heavier singleton, with intersection volume
   (L/2-7/16)*L^4. The exact answer is

       37/64 + (7/16)*(31/32)^(4/5)
       = 1.004652885572959045694... .

The full rational enclosures, execution timings, node counts and hashes are
in the three result files and `completion_receipt.json`. No floating-point
sample is used to decide a bound. The symbolic residual identities are also
actually checked with the installed SymPy runtime in `proof_identity_check.py`.

## Reproduction and precise next question

Run the three commands `python3 two_atom_cluster.py equal_axis`,
`python3 two_atom_cluster.py equal_diagonal` and
`python3 two_atom_cluster.py unequal_axis` from this directory. Each input
has the independent 120-second / 200,000-node limit. Run the identity checker
with the existing local Python environment containing SymPy, whose actual
path and version are recorded in its output. No MATLAB is used.

The unresolved next finite problem is to choose three distinct points and
positive weights, list all seven subset masses and common-radius receiver
boxes, and prove/compute their union ratio. This round does not launch that
search, claim a promising triple, or suggest that testing only pairs can
establish a general C5 upper bound.
