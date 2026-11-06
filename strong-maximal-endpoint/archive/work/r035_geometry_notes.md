# r035 actual cross-rectangle finite-grid probe

Status: deterministic finite-grid evidence only.  It does not certify a
continuous two-dimensional optimum, improve (C_2), or lift to correlated
high dimension.

The input is constant on the sixteen unit cells of ([0,4]^2).  Receiver
cells use edges `[-4,-1,0,1,2,3,4,5,8]`.  All 1296 rectangles whose four
edges lie in this list are enumerated.  A receiver cell is counted only when
it lies wholly in one rectangle with exact integer test

\[
1000\,(\text{input mass})>1001\,(\text{rectangle area}).
\]

The final receiver set is a bitwise union over these actual valid rectangles.
It is not a sum of component halos.

The search covers every four-cell support and amplitudes (2,4,8,16).  Its
best disconnected ratio occurs at amplitude (4), with support

\[
\{(1,2),(1,3)\}\cup\{(3,2),(3,3)\}.
\]

There are two 1-by-2 disconnected components.  The union from the full input
has grid area 48; a **finite-template diagnostic** obtained by separately
evaluating the components with the same 1296 rectangles has area 38.  Its
difference 10 records a cross-witness effect inside this finite template, but
is not the difference from the full true component halos.  The exact Orlicz cost is

\[
4\cdot4(1+\log4)=16+32\log2,
\]

and the grid score is (48/(16+32\log2)=1.2571793525892\ldots).
On the same receiver grid, the uniform (4,1_{[0,4]^2}) template has union
area 108, cost (64+128\log2), and score (0.7071633858314\ldots).
This is a same-grid comparison, not a continuum comparison.

The six cells responsible for the additional area, and a specific rectangle
for each, are in `r035_geometry_witness.json`.  For example,
`[-4,3] x [1,3]` captures all four high cells and has exact average (8/7),
strictly exceeding (1001/1000).

The template probe is outside the scope of the r007 component-wise union
calculation in the limited sense that it admits rectangles taking mass from
more than one component.  It nevertheless has no automatic high-dimensional consequence.
The separate r035 tensor-product obstruction shows that naively repeating a
bounded two-dimensional base incurs a denominator cost of order
(\Theta(nk\epsilon)).  A possible next model would need a genuinely
non-product block construction with cross-block rectangles, together with a
dimension-uniform accounting of the true Orlicz cost and receiver union.

Separately, `r035_geometry_certificate.json` is a genuine continuous,
all-rectangle statement.  For the four isolated cells
`(1,2),(2,3),(3,2),(4,3)` at amplitude 16, the square

\[
E=[-4,-7/2]^2
\]

has area (1/4).  The rectangle ([x,4]\times[y,3]) captures all mass 64
for each ((x,y)\in E), with average at least (8/7).  In contrast, the
exact strong-maximal values of the four individual unit cells are at most
(64/99,64/143,64/143,64/195) on (E), all below one.  Thus (E) is in the
full input halo but no individual component halo.  The binary (4\times4)
support matrix has rank 2, so this particular cross witness is non-product.
