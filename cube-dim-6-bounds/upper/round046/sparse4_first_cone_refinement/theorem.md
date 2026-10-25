# A046: complete exact refinement of sparse4 parent cone1

## Scope and normalization

The new frame is exactly the B044 locked 152-column matrix A=[body32,(1/80)W], W the 120 unoriented support-four sign columns. Its characteristic exponent is

    ψ(ξ)=Σ_(a∈A)|a·ξ|.

All columns, including the factor1/80 on W, are read from the locked frame.json. This theorem proves only the refinement of old parent cone1. It does not construct the other20 parent refinements, a global fan, a heat kernel, or a Cn bound. No old tube or value is transferred.

Let R have columns

    r1=(2,1,1,1,1,0), r2=(3,3,2,2,2,0), r3=(2,2,2,1,1,0),
    r4=(1,1,1,1,0,0), r5=(1,1,1,1,1,0), r6=(1,1,1,1,1,1).

Exact detR=1. Every nonzero point of this cone is a positive multiple of Ru for a unique u≥0 with Σu=1. Thus a complete subdivision of the standard5-simplex P={u≥0,Σu=1} gives a complete cone subdivision, with the origin included by closure. Normalized simplex volume means volume relative to P, whose u-coordinate Euclidean5-volume is √6/120.

The physical slice has equation ℓ·ξ=1 with ℓ=R^−T1=(0,−1,1,1,0,0), |ℓ|²=3. Its Euclidean5-volume is √3/120, not √6/120. The ratio of sub-simplex to parent volumes is the same in both coordinate systems. For six slice vertices u^1,…,u^6,

    relative volume = |det[u^1 … u^6]|.

This determinant formula follows by replacing one affine coordinate with the common sum coordinate. It does not use a floating hull or a numerical volume.

## 1. Independent first-cone crossing check

Compute only the120 new forms on these six rays, 720 integer dot products. The crossing indices, in the locked W order and using1-based indexing, are exactly

    2,10,26,50,82.

Their pullbacks to u coordinates are

    h1=u1−u3,
    h2=u1−u3−u4,
    h3=u1−u4,
    h4=u1+u2−u4,
    h5=u2−u4.

Each has both strict signs at parent vertices. The other115 new forms have one weak sign on the parent, including zero contacts. The old2160-ray/form baseline was not rerun. These five cuts all ignore u5,u6. In particular Pσ below is the join of its first-four-coordinate slice with the edge conv(e5,e6). This explains why the nontrivial refinement is three-dimensional despite the5D slice; no such simplification is assumed for the other20 parents.

## 2. Exact halfspace partition and complete vertices

For every σ∈{−1,+1}^5 define the closed cell

    Pσ={u:Σu=1, ui≥0 (i=1,…,6), σj hj(u)≥0 (j=1,…,5)}.

There are32 registered patterns and11 halfspaces per cell. Their union is P: choose each sign to match hj(u), choosing either at a zero. If two patterns differ, their intersection lies in at least one of the cut hyperplanes. Thus distinct full cells have disjoint relative interiors; zero/tie faces are retained.

Every vertex of a bounded polyhedron in the affine5-space Σu=1 has at least5 independent active boundary normals within that affine space. Hence it occurs among the systems consisting of Σu=1 and five of the11 boundary equalities. Conversely a feasible solution to an independent such system is a vertex. Boundary equalities do not depend on σ, so the matrices can be solved once and their feasible points filtered separately for every pattern. This is a cache of exact systems, not an assumption about floating candidates.

The actual enumeration exhausts binom(11,5)=462 matrices:237 are singular,6 nonsingular solutions are outside u≥0, and the remaining solutions deduplicate to12 global vertices. Each saved vertex has an independent active-set witness. Feasibility against all5 signed cuts is tested exactly for all32 patterns. The result is12 cells of affine dimension5 and20 nonempty lower-dimensional cells; no pattern is silently omitted. They all contain the common edge conv(e5,e6).

Each lower-dimensional pattern's saved complete vertex set lies in one recorded full-dimensional closed pattern. Taking convex hulls proves its entire cell lies there. More generally, full-dimensional arrangement cells cover P by closure: the complement of finitely many nonzero cut hyperplanes is dense in P, each such point lies in a full cell, and a finite union of closed cells is closed. Either argument retains all boundary, tie and zero cases.

## 3. Exact facets and recursive pulling

For a recorded polytope face F with complete vertex set S and affine dimension d, every facet is obtained by imposing equality in at least one original defining halfspace. Intersect S with each of the11 boundary hyperplanes; retain precisely distinct proper intersections of affine dimension d−1. This finds all facets, including nonsimplicial or multiply supported facets. All affine dimensions are exact rational ranks.

Use one global lexicographic order of all12 rational u vertices. The pulling subdivision of F is defined inductively. If d=0 it is the vertex itself. Otherwise let p be its least vertex. For each facet E not containing p, recursively triangulate E and cone every resulting simplex to p.

Here is the coverage and non-overlap justification. From p, a ray through an interior point of F exits F through a facet not containing p, except on lower-dimensional ties. It follows that F is the union of the pyramids conv(p,E) over those facets. Their interiors are disjoint because the first exit facet is unique away from facet intersections. Each pyramid is subdivided by the inductively complete triangulation of its base. Closure retains every boundary point. This proves coverage and disjoint interiors by dimension induction, rather than inferring coverage from a floating hull or merely trusting the final volume.

The least-vertex rule restricts consistently to faces, so two incident cells obtain the same triangulation on their intersection. In addition to this induction, the construction checks all66 pairwise nonempty intersections of the12 full cells: the maximal induced simplices on each shared face agree exactly. All66 recursively used face nodes, their facets, apex and simplices are saved as a proof DAG in local_fan.json. A second checker reconstructs these facet/rank relations and every recursion identity from the saved certificate, without re-enumerating active systems.

## 4. True primitive rays, linear rates and determinant weights

For a rational slice vertex u, set x=Ru and uniquely write x=s r, s>0 rational, r primitive integer. The output records u, the active witness, r and s. Primitive rays remain in the ordered Fourier chamber because the parent rays do and u≥0.

For each r the rate is recomputed from all152 true columns,

    q(r)=ψ(r)=(1/80)Σ_(j=1)^152 |column_numerator_j·r|.

No old q is read or reused. Every full cell's complete vertex list has one weak sign for each of the152 forms; linearity then gives the same sign throughout the cell. If εj is that sign, the exact slope is

    bσ=(1/80)Σ εj column_numerator_j,
    ψ(ξ)=bσ·ξ   on cone(Pσ).

The saved linear coefficient is independently checked against q(r) at every child ray. Every output simplex also passes the form-sign test separately on all its six primitive rays. Positive scaling extends this sign and linearity assertion to its entire closed cone.

For each six-ray child, exact det[u^1 … u^6] is nonzero. If x_i=s_i r_i, then

    |det[r_1 … r_6]| ∏s_i=|detR| |det[u^1 … u^6]|.

Both determinant sides are checked, giving the correct primitive physical determinant, not a determinant of untracked normalized rays.

There are14 output full-dimensional simplices. Their normalized volumes are

    1/6,1/24,1/8,1/24,1/8,1/12,1/48,
    1/48,1/24,1/48,1/48,1/24,1/12,1/6,

with exact sum1. Thus the u-coordinate volumes sum to √6/120 and physical slice volumes sum to √3/120. Volume equality is a redundant exact check of the partition/pulling proof, not its replacement.

## 5. Acceptance object and limits

local_fan.json contains all32 closed patterns,12 primitive rays with true rates/scales,66 face proof nodes, and14 simplicial cones with exact physical determinants and volume ratios. This is a complete legal local fan for the first old parent cone at beta1/80. Its shared faces within this parent are compatible. It makes no claim about a full frame fan or the yet-unrefined neighboring parents.

The rest of a valid global Fourier interface still requires the other20 parent refinements and their boundary integration, global symmetry multiplicities/normalization, and treatment of repeated/zero Fourier nodes. No old21-cone q-only engine is justified by this one-parent result. No heat-kernel positivity/numerical value, time budget, singular tube or Cn improvement is certified here.
