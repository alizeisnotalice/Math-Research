# Exact two-orbit witness transfer

This proves a lower bound from either rational witness file; neither exceeds
the standing C. It does not identify the full periodic superlevel set.

For the recorded rational spacing vector a, put unit mass at diag(a)k and
mass ε=1/1024 at diag(a)(k+b), k∈Z⁵, b=(0,1/2,0,0,0). The two orbits are
disjoint: their second coordinates have distinct integer/half-integer parity.
Every period has mass 1+ε. No atom is counted twice in a selected witness.

For each recorded witness W, its exact mass m_W is the sum of these distinct
atomic masses and its common radius r_W is rational with (2r_W)⁵<m_W.
The receiver box is the intersection, for all five coordinates, of

    [max_{p∈W} p_j-r_W, min_{p∈W} p_j+r_W]

with H=∏[0,a_j/2]. The two endpoint inequalities have been checked against
every selected atom. Thus every point of that box is the centre of one
five-dimensional cube, with the same radius r_W in all coordinates,
containing the entire selected mass and having average strictly above one.
The exact union volume V is counted once, including all overlaps.

Both orbit measures are invariant under reflection in each coordinate and
translation by diag(a)Z⁵. Reflecting H across its midpoint faces covers
∏[0,a_j], up to null boundaries. Hence one full period receives volume 32V,
and the finite-witness ratio is A=32V/(1+ε). There is no further division by
the geometrical period volume.

The original witnesses use lattice indices of absolute value at most two.
Reflection x_j↦a_j-x_j sends a primary index k_j to1-k_j and a satellite
index to1-k_j-b_j-b_j; this is integral because2b_j∈{0,1}.
Therefore all reflected witnesses have index coordinates between-3 and3.
For M≥1 retain both orbits with indices k∈{-3,…,M+3}⁵, and receive only in
the M⁵ full periods whose indices are in{0,…,M-1}⁵. This collar contains every
needed atom. The exact input mass is (1+ε)(M+7)⁵, and the received volume is
32V M⁵. The finite atomic weak ratio is at least

    A [M/(M+7)]⁵.

To obtain genuine smooth L¹ inputs, use the strict radius inequalities.
There are finitely many witnesses, so choose δ>0 so small that
(2(r_W+δ))⁵<m_W for every W. Replace each retained atom by a nonnegative
smooth bump of its same mass supported within a cube of radius δ/2 around
that atom. Every recorded receiver cube enlarged to common radius r_W+δ
contains the entire selected bumps and still has average strictly above one.
The mass is unchanged; the bumps have compact support and the sum is smooth
and nonnegative. Thus the same receiver-volume bound is valid for L¹ inputs.
Letting M→∞ yields C₅≥A at threshold one. The maximal operator includes all
positive common radii, so the finite witness scales are admissible. This
argument uses fixed Lebesgue background throughout.
