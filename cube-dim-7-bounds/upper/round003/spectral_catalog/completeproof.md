# Round 003: complete seven-dimensional spectral vertex orbit catalog

## Object and result

Fix the **unchanged** seven-dimensional spectral variable `xi` and the 64 body
directions

\[
 V=\{(1,s_1,\ldots,s_6):s_j\in\{-1,1\}\},\qquad
 \psi(\xi)=\sum_{v\in V}|v\cdot\xi|,\qquad
 P=\{\xi:\psi(\xi)\le1\}.
\]

This is the original same-time spectral model. There is no change of time
coordinate and no passage to physical Green variables.

**Certified result.** The vertices of `P` have exactly **21 orbits** under all
signed coordinate permutations, with **504868 vertices** in total. Every orbit
appears in `catalog.json` with its primitive integral direction, exact integral
`psi_rate`, all seven reduced rational coordinates of the unit-ball vertex,
all zero-direction indices, an explicit independent six-row witness, a nonzero
minor, all seven signed cofactors, and the full orbit cardinality.

This proves a finite spectral classification. It proves no assertion about
the physical Green singular set, physical `S3`, a minimum, `g`, or `C_n`.

## 1. Nondegeneracy and the exact vertex criterion

The 64-by-7 matrix of rows `V` satisfies `V^T V=64 I`: the product of two
different coordinate columns sums to zero over the six independent signs,
and every column has squared length 64. Thus `V` spans `R^7`, and `psi` is a
norm. In particular `P` is a full-dimensional bounded polytope. Indeed

\[
 \psi(x)=\max_{\varepsilon\in\{-1,1\}^{64}}
       \left(\sum_{v\in V}\varepsilon_v v\right)\cdot x.
\]

The last representation is used only as a mathematical identity; no such
sign-family enumeration is run.

For nonzero `x` write

\[
 Z(x)=\{v\in V:v\cdot x=0\},\qquad
 b(x)=\sum_{v\cdot x\ne0}\operatorname{sgn}(v\cdot x)v.
\]

The supporting gradients active at `x/psi(x)` are exactly

\[
 A(x)=\left\{b(x)+\sum_{v\in Z(x)}\varepsilon_v v:
       \varepsilon_v\in\{-1,1\}\right\}.
\]

The signs at a nonzero dot product are forced by equality in the absolute
value identity; signs at the zero dots are free. Flipping a single free sign
shows that

\[
 \operatorname{aff} A(x)=b(x)+\operatorname{span}Z(x),\qquad
 \operatorname{span} A(x)=\operatorname{span}(b(x),Z(x)).
\]

Since every `v` in `Z(x)` is perpendicular to `x` but
`b(x) dot x = psi(x)>0`, the latter span has dimension `rank Z(x)+1`.
Also `rank Z(x)<=6`.

If `rank Z(x)=6`, then the active gradients span `R^7`. Suppose the
unit-ball point is the midpoint of `y,z in P`. Every active supporting
gradient has value 1 on the midpoint and at most 1 on either endpoint, so
has value 1 on both endpoints. The active gradients spanning `R^7` force
`y=z`. Hence the point is a vertex.

Conversely, if `rank Z(x)<6`, choose a nonzero vector `h` perpendicular to
`b(x)` and to every member of `Z(x)`. Such a vector exists because their
combined span has dimension at most 6. Normalize `x` to `psi(x)=1`. For a
sufficiently small positive `t`, all nonzero `v dot x` keep their sign in
`v dot (x +/- t h)`, while zero dot products remain zero. Therefore

\[
 \psi(x\pm t h)=b(x)\cdot(x\pm th)=1.
\]

These are distinct boundary points with midpoint `x`, so `x` is not a
vertex. This supplies the asserted two-sided boundary line segment and
proves the equivalence

\[
 \boxed{x/\psi(x)\text{ is a vertex of }P
       \iff\operatorname{rank}Z(x)=6.}
\]

## 2. Integer cofactor bound and completeness of the finite domain

For any vertex direction select six independent rows of `Z(x)` and call
the resulting 6-by-7 matrix `B`. Its entries are all `+/-1` and its rank
is 6. The signed 6-by-6 maximal minors

\[
 c_j=(-1)^j\det B_{\widehat j},\qquad j=0,\ldots,6,
\]

form a nonzero vector in the one-dimensional kernel of `B`. Thus `x` is
proportional to `c`.

For **each** 6-by-6 matrix with entries `+/-1`, subtract the first row
from each of the other five rows. Every entry in each of those five rows
is even; multilinearity of the determinant gives divisibility by
`2^5=32`, including when the determinant is zero. Hadamard's determinant
inequality gives

\[
 |\det B_{\widehat j}|\le(\sqrt6)^6=216.
\]

Hence `d=c/32` is a nonzero integral kernel vector and every coordinate
is an integer of absolute value at most `floor(216/32)=6`. In fact its
unscaled maximal minor cannot have absolute value between 193 and 216
because it is a multiple of 32; no sharper determinant theorem is needed.

Divide `d` by the positive gcd of its coordinates. This produces a
primitive integral direction `p`; taking the gcd cannot increase any
coordinate. Therefore every vertex direction has a primitive integral
representative with all coordinate absolute values in `{0,...,6}`.

This argument proves completeness before computation. It needs no scan
over six-element subsets of the 64 body directions, no old baseline, and
no hyperplane fan computation.

## 3. Symmetry, unique orbit representatives and orbit sizes

Let `Q` be any signed permutation matrix. Applying `Q^T` to `V` sends
each body direction to a sign times a unique member of `V`: the full
sign cube consists of the 64 represented antipodal pairs, and `Q^T`
permutes these pairs bijectively. Thus

\[
 \psi(Qx)=\sum_{v\in V}|(Q^Tv)\cdot x|=\psi(x).
\]

The action preserves `P`, its vertices, and the rank criterion. Signs and
permutations send every primitive integral direction to a unique
nonnegative nonincreasing tuple of seven integers with gcd 1. Conversely,
two such canonical tuples on the same positive ray are equal: primitive
integral points on a rational ray are unique. Equal unit-ball vertices
imply equal rays, so no two different canonical tuples describe the same
orbit of vertices.

The number of nonincreasing seven-tuples drawn from the seven-element
alphabet `{0,...,6}` is the multiset count

\[
 \binom{7+7-1}{7}=\binom{13}{7}=1716.
\]

The zero tuple and nonprimitive tuples are excluded. For each remaining
tuple the exact rank criterion is both necessary and sufficient, so this
finite procedure enumerates **all and only** the vertex orbits.

For a representative `p`, let `k` be its number of nonzero coordinates,
and `m_a` the multiplicity of each distinct coordinate value, including
zero. Its distinct signed permutations number

\[
 |\mathcal O(p)|=2^k\,\frac{7!}{\prod_a m_a!}.
\]

There are `7!/prod m_a!` distinct coordinate placements and, independently,
`2^k` choices of signs on the occupied nonzero positions. This includes
the antipodal vertex; zero signs create no additional vertices. Since
`psi` is invariant, normalizing by `psi(p)` does not alter this count.
Different canonical representatives have disjoint orbits, so summing
these integers gives the exact total vertex count.

## 4. Exact computed catalog

The following is a compact view of `catalog.json`. `z` is the number of
zero body directions; every displayed row has exact zero-direction rank 6.
Coordinates of the corresponding unit-ball vertex are exactly `p/rate`.

| ID | primitive p | rate | z | orbit size |
|---|---|---:|---:|---:|
| 001 | (5,4,3,2,2,1,1) | 404 | 6 | 161280 |
| 002 | (5,3,3,2,1,1,1) | 372 | 6 | 53760 |
| 003 | (5,2,2,2,1,1,1) | 340 | 6 | 17920 |
| 004 | (4,3,3,2,2,1,1) | 344 | 7 | 80640 |
| 005 | (4,3,3,1,1,1,1) | 324 | 6 | 13440 |
| 006 | (4,3,2,2,1,1,1) | 312 | 8 | 53760 |
| 007 | (4,2,2,1,1,1,1) | 280 | 8 | 13440 |
| 008 | (4,1,1,1,1,1,1) | 260 | 6 | 896 |
| 009 | (3,3,2,2,2,1,1) | 292 | 8 | 26880 |
| 010 | (3,3,2,1,1,1,1) | 264 | 9 | 13440 |
| 011 | (3,2,2,2,1,1,1) | 252 | 10 | 17920 |
| 012 | (3,2,2,1,1,1,0) | 232 | 10 | 26880 |
| 013 | (3,2,1,1,1,1,1) | 220 | 11 | 5376 |
| 014 | (3,1,1,1,1,1,0) | 200 | 10 | 2688 |
| 015 | (2,2,2,1,1,1,1) | 204 | 12 | 4480 |
| 016 | (2,2,1,1,1,1,0) | 176 | 14 | 6720 |
| 017 | (2,1,1,1,1,1,1) | 160 | 15 | 896 |
| 018 | (2,1,1,1,1,0,0) | 144 | 16 | 3360 |
| 019 | (1,1,1,1,1,1,0) | 120 | 20 | 448 |
| 020 | (1,1,1,1,0,0,0) | 96 | 24 | 560 |
| 021 | (1,1,0,0,0,0,0) | 64 | 32 | 84 |

The 1716 tuples split into 1 zero tuple, 154 nonzero nonprimitive tuples,
and 1561 primitive tuples. The latter have zero-direction rank histogram

| rank | 0 | 1 | 2 | 3 | 4 | 5 | 6 |
|---|---:|---:|---:|---:|---:|---:|---:|
| count | 954 | 13 | 30 | 109 | 242 | 192 | 21 |

Thus all 1540 rejected primitive tuples have exact rank at most 5. The
code regenerates each of them and its rank; storing the histogram is
sufficient and avoids unnecessary repeated data. The complete domain
included coordinate value 6; none of the accepted primitive tuples has
maximum coordinate 6. This observation is a computed result, not a
replacement for the proved bound 6.

## 5. Computation, receipt and independent replay

`registration.json` was written before the first mathematical run. The
registered work bound is the finite domain of at most 1716 tuples, not
a wall-time cutoff. `catalog.py` performs 99904 integer dot products on
the 1561 primitive candidates and calls `python-flint fmpz_mat.rank` on
their exact zero-direction matrices. For accepted cases only, exact
SymPy RREF chooses six independent rows and six pivot columns, then
exact determinants provide the minor and cofactor witnesses. All
coordinates are serialized as reduced rational strings, never decimals.

The actual primary run used Python 3.15.0, python-flint 0.9.0, and
SymPy 1.14.0. It took approximately 0.1023 seconds wall time and 0.1015
seconds process CPU. These measurements are recorded in `receipt.json`;
no multicore utilization claim is made. The receipt records source and
result SHA-256 hashes. `run.stdout` is the actual short successful stdout,
and `run.stderr` is empty.

`recheck.py` is a **stdout-only** checker: it never writes a file, uses
no third-party package, and calls no floating-point or library rank
routine. It independently reconstructs all 1716 tuples, evaluates zero
forms, and determines exact rank using integer row elimination with gcd
reduction. Each operation replaces a lower row `R_i` by
`a_rc R_i - a_ic R_r`; since the nonzero pivot `a_rc` is invertible over
the rationals, this preserves the rational row space. Dividing a nonzero
row by its gcd also preserves it. Pivot counting therefore yields the
exact rank. It checks that no generated integer exceeds the registered
768-bit ceiling. Its six-by-six determinant checker uses the independent
Leibniz formula, verifies every stored minor and cofactor, every exact
unit-ball coordinate and norm, every orbit size, all histogram counts,
and the total. The measured maximum integer width for its elimination
intermediates is reported with its actual stdout; no hidden library
allocation width is claimed. The independent replay does not perform
any baseline or Green calculation.

Run from any directory, with the existing environment:

```text
/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python -B /Users/zhengzhihao/Desktop/T/Math-Research-cube7-spectral-r003/cube-dim-7-bounds/upper/round003/spectral_catalog/catalog.py
/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python -B /Users/zhengzhihao/Desktop/T/Math-Research-cube7-spectral-r003/cube-dim-7-bounds/upper/round003/spectral_catalog/recheck.py
```

The first command recreates the primary outputs and receipt in this
directory. The second only reads files and prints its result. The saved
`recheck_receipt.json` records a real subprocess invocation of the second
command and its stdout. `manifest.json` lists the final file whitelist,
sizes and hashes; its own hash is intentionally not embedded in itself.

## 6. Provenance and proof status

The standalone clone was created with `--no-hardlinks --no-checkout`
from the local source clone at
`b2139fc1fe70fd55e7d409031256195c7d10e768`, then configured with the
requested SSH origin, `promisor=true`, `partialclonefilter=blob:none`,
`extensions.partialClone=origin`, and sparse path `cube-dim-7-bounds`.
Only the present `upper/round003/spectral_catalog/` directory is written.
No commit, push, remote fetch, remote acceptance claim, or outside-scope
cleanup was performed. `sources.json` records the local source identity.

Sections 1--3 are the independent mathematical audit of the finite
reduction supplied in the task. They are self-contained proofs, using
elementary linear algebra, the determinant row operation, Hadamard's
inequality, and multiset counting. Sections 4--5 provide an exact finite
computation and independently replayed certificates. There are no
remaining gaps for this spectral catalog statement. Transfer from this
catalog to physical Green geometry or a new upper bound remains outside
this result and would require a separate argument.
