# A012: exact five-term base and six-term augmented distance formulas

Both proposed distance identities are **true**. The augmented eight-term
version has two redundant terms and simplifies to **six**. The proof is
global on R6, including signed coordinates, ties and zeros, and does not
rely on finite enumeration. The finite rational verification then checked
all496 base pairs and all703 augmented pairs at each of twelve points.

For z=sort_desc(|x|), put

    W(a1,...,am)=sum a_i^2-(sum a_i)^2/m,
    E_k=W(z1,...,zk)+W(z_(k+1),...,z6),  k=1,...,5.

The exact formulas are

    dist(x,S1(body32))^2=W(z1,...,z6),
    dist(x,S2(body32))^2=min(E_1,...,E_5),
    dist(x,S2(body32+axes))^2
       =min(E_1,...,E_5,z3^2+z4^2+z5^2+z6^2).

For completeness, the augmented nearest-line distance is
min(W(z1,...,z6),z2^2+...+z6^2). Every axis coefficient beta>0 has the same
span geometry. These formulas are geometric statements, not G/K or operator
certificates.

## Why five contiguous partitions suffice, including all ties

The full signed-permutation group preserves the32 unoriented body lines,
their pair spans and Euclidean distance. Hence reduce to nonnegative sorted
z. Any body-pair plane splits its coordinates into two nonempty blocks on
which the two body directions agree or disagree. Its orthogonal basis is
the two disjoint signed block indicators. The squared projection residual is

    ||z||^2-(signed block sum_A)^2/|A|
             -(signed block sum_B)^2/|B|.

Absolute signed sums are at most the corresponding unsigned sums. Every
unsigned block partition is realized by the all-positive body direction
and a body direction that changes sign on one block. Thus the exact nearest
plane problem is the nonempty two-cluster within-variance problem.

A fixed-cardinality argument removes any nearest-center reassignment
ambiguity. For |A|=k, S=sum z_i and s=sum_A z_i, maximizing projection gain
means maximizing

    f_k(s)=s^2/k+(S-s)^2/(6-k).

This is convex, while s lies between the bottom-k sum and top-k sum. Its
value is at most the greater endpoint value. Thus one may take the top k
or bottom k; the latter is the complement of top6-k. The five top splits
therefore suffice. This proof never empties a cluster, and ties do not
change endpoint sums. Equal/zero values and nonunique closest planes are
all included. The full self-contained proof is in `theorem.md`.

For implementation, if S_k=sum_(i<=k)z_i and Q6=sum z_i^2, each candidate is

    E_k=Q6-S_k^2/k-(S6-S_k)^2/(6-k).

Only sorting six absolute values and five residual evaluations are needed.

## The augmented leave-one-out terms are already base terms

The best body-axis plane has one free coordinate and a constant signed
line on the other five. Removing value t gives

    L(t)=sum z_i^2-t^2-(S-t)^2/5,
    L''=-12/5<0.

Its minimum over the observed values is attained by deleting a maximum
or minimum. The two proposed terms are therefore W(z2,...,z6) and
W(z1,...,z5). They equal E_1 and E_5 respectively.

Geometrically, if w is v with coordinate i flipped,

    w=v-2v_i e_i,
    span(v,e_i)=span(v,w).

Thus body-axis planes add nothing to base S2. Axis-axis planes add the
single nearest-distance candidate obtained by retaining the top two absolute
coordinates. This proves the eight-term proposal and its six-term reduction
for every beta>0.

The exact full-frame projector comparison also found:

|Pair class|Direction pairs|Geometrically distinct planes|
|---|---:|---:|
|Body-body|496|496|
|Body-axis|192|96, all already among the base496|
|Axis-axis|15|15, all additional|
|Augmented total|703|511|

The distinct-plane count has an analytic proof as well: each base pair plane
contains exactly its two body lines, and every coordinate2-plane contains
no body direction. The union is the old496 planes plus15 coordinate planes.
This union redundancy does not remove axes or their weighted pair terms
from a density/product estimate; beta still affects the kernel and its
analytic constants.

## B010 exact inputs reproduced

The two read-only B010 rational inputs were used directly, rather than their
printed decimal coordinates. Their original exact squared S2 distances were
recovered by both the five-term formula and exhaustive physical projection:

    point1: 21826641225836976227261
              /200000000000000000000000,

    point2: 54771921054251644984225195390456050322661189
              /500000007000000024500000000000000000000000000.

Their distances have Arb384 enclosures around

    0.330353153654062883590075412887,
    0.330974078403263352228037084566.

Both attain the contiguous split k=4, and their augmented S2 distances remain
the same. The physical closest base pair indices28,31 also match the B010
enumeration. No G/H/K value was evaluated or transferred to another frame.

## Finite execution and useful counterchecks

Exactly twelve rational vectors were registered before the complete-pair
checks. They include origin, equal and signed equal values, a single axis,
two unequal axes, repeated two/three levels, a signed/permuted control, a
strict chamber control, a signed body-axis vector and the two B010 inputs.

Every point matched all496 base pair projections, all703 augmented pair
projections and the corresponding32/38 nearest-line projections with exact
Fraction arithmetic. The augmented test used the fixed rational beta=3/7;
the analytic proof handles every positive beta. No second beta or extra
point was tested.

For instance, x=(1,1/3,0,0,0,0) has base squared distance4/45 and augmented
distance0. This confirms that the extra coordinate planes genuinely enlarge
S2, although the body-axis planes do not. The repeated-level point
(1,1,1/2,1/2,0,0) has two optimal contiguous splits k=2 and4; the exact
minimum1/4 is correctly retained.

The one script completed with exit0 in **0.148197917 seconds** (tool chunk
bfbc72), Python PID **28499**. SymPy additionally checked
the convex gain and concave leave-one-out derivatives, and the dependent
augmented triple. There were zero kernel calls, MC samples, MATLAB calls,
agents, commits or pushes. All input definitions and source hashes are in
`input_manifest.json`, and every exact candidate/minimum is in `results.json`.

## Coverage and rank cautions

The five base formulas calculate the nearest distance in a sorted absolute
chamber. They are not a replacement of global S2 by five fixed physical
planes, and no pair-plane tubes are asserted disjoint. At ties many more
physical planes can be nearest; for example the all-equal body point belongs
to31 base pair planes. A deterministic tie-break can select an argmin, but
does not justify adding tube volumes as if they were a disjoint partition.

Three distinct augmented directions can be dependent: with
v=(1,...,1), w=(1,-1,1,1,1,1) and beta e2,

    beta v-beta w-2(beta e2)=0.

The symbolic rank is exactly2. Future large-variable decompositions must
use rank-greedy independent subsets; the selected basis spans the full large
direction set and its size is the actual rank. Cardinality alone cannot
replace rank. No new-frame singularity classification or covering constant
is inferred here.

This finite task is complete and stopped. The six retained files are the
registration, theorem, this report, checker, exact results and input manifest.
All previously completed directories are unchanged. Reproduce with:

    /Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python check_geometry.py
