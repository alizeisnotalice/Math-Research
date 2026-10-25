# Exact projection peaks by a one-sided full marginal trace

The original body32 symbol is psi(xi)=sum_v |v.xi|. Its complete fan is
read from the existing exact21-cone certificate, together with all720
coordinate permutations and64 signs. Cone interiors partition R6 up to
Lebesgue null sets. The original exact coverage proof is reused, not
replaced by a numerical fan or by boundary counting. In every base cone,
xi=R u,u>=0, psi(xi)=q.u, and |det R|>0. Its normalized weight is
alpha=|det R|/product_j q_j.

Let L have r=1,2,3 independent selected body directions as rows, W=ker L,
d=6-r, and D=det(LL^T)>0. Let q_d be the density of the orthogonal projection
of X_1 onto W. The selected directions themselves vanish in this quotient,
so this is also the quotient law of the remaining directions in B008.
With orthonormal Hausdorff/Lebesgue measure on W,

    q_d(1,0)=(2pi)^-d integral_W exp(-psi(w)) dH^d(w).            (1)

The symbol restricted to W is coercive and positive homogeneous; hence
the integral is finite. Fourier inversion and the triangle inequality give
q_d(1,y)<=q_d(1,0); positivity follows from the actual projected Cauchy law.
Scaling gives sup_y q_d(t,y)=t^-d q_d(1,0) for every t>0.

## 1. Use the full marginal before restricting to a fan boundary

Push the measure exp(-psi(xi)) dxi forward by L. Its density is

    F(z)=D^-1/2 integral_W exp[-psi(w+L^T(LL^T)^-1 z)] dH^d(w). (2)

The coarea Jacobian is sqrt(D); it is not its reciprocal in(1).
By the triangle inequality for psi and coercivity, on bounded z sets
this integrand is bounded by a constant times exp(-c|w|). Dominated
convergence proves that F is continuous everywhere, including at0.
Thus integral_W exp(-psi)=sqrt(D)F(0).

The cone decomposition of the original six-dimensional measure pushes
forward as an equality of measures. We will evaluate the cone density
versions at a point where all are continuous, then pass to0. We do NOT
partition W directly by the fan: W lies on body-form walls and such a
partition could double-count faces of full intrinsic measure.

## 2. One cone and the number of nonzero normal columns

In one cone put y_j=q_j u_j. The measure becomes
alpha exp(-sum y_j)1_(y>=0)dy. Its normal map is

    A=L R diag(1/q_j).

Let J be the nonzero columns, m=|J|. Since rank L=r and R is invertible,
rank A=r and m>=r. The zero-column variables integrate out to1. Every
selected body form has constant weak sign on the cone, so each row of A
has no mixture of positive and negative entries. After individually
flipping rows, all columns are nonnegative; every column in J is nonzero.
Consequently there is a row sign vector s such that s.a_j>0 for every
j in J. The nonnegative fiber A_J y=z is bounded: if it is nonempty,

    sum y_j <= (s.z)/min_j(s.a_j).

This prevents an unbounded fiber from cancelling the positive scaling
power used below.

If m=r, M=L R_J is invertible, and the cone pushforward at z is

    alpha/|det(A_J)| exp[-sum(A_J^-1 z)] 1_(A_J^-1 z>0)         (3)

away from its boundary. Since q_j>0, the sign of A_J^-1 v is exactly the
sign of M^-1 v. At z=epsilon v, epsilon down to0, the limit is zero if
M^-1v has a negative component, and otherwise is

    alpha/|det(A_J)|
      =|det R|/[|det M| product_(j notin J)q_j].                (4)

If m>r, coarea on A_J and y=epsilon w give the density along the ray as

    alpha epsilon^(m-r)/sqrt(det(A_J A_J^T))
      * integral_{A_Jw=v,w>=0} exp(-epsilon sum w_j) dH^(m-r)(w).

The bounded fiber makes the integral uniformly finite. Hence this
contribution tends to0. If the fiber is empty, the density is zero.
There are finitely many cones, so no infinite limit interchange occurs.

## 3. Why equality of densities holds on the selected ray

Choose v outside every hyperplane spanned by r-1 independent nonzero
normal columns of every cone. For r=1 require v!=0. This also excludes
m=r coordinates of M^-1v being zero.

For completeness, each cone pushforward is continuous near epsilon v.
Choose any independent r columns B. Eliminate their coordinates, leaving
m-r free nonnegative coordinates w; the density is an ordinary integral
with factor1/|det A_B| and indicators
A_B^-1(z-A_free w)>=0. The bounded-fiber estimate gives a uniform compact
integration domain near z=epsilon v. A nonconstant affine boundary in w
has zero (m-r)-dimensional measure. If a solved coordinate is independent
of w, its vanishing is a hyperplane in z spanned by r-1 columns, excluded
by the choice of v. Thus dominated convergence gives continuity here.
The m=r case is directly continuous away from the same hyperplanes.

Both sides of the pushed-forward measure identity are therefore continuous
in a neighborhood of epsilon v; their almost-everywhere equality extends
to this point. Now F is continuous at0. Summing(4) yields the exact number

    R_L := F(0)
       =sum_{cones:m=r, M^-1v>0}
          |det R|/[|det M| product_(j notin J)q_j].              (5)

Finally,

    q_d(1,0)=sqrt(D) R_L/(2pi)^d,
    kappa_L=q_d(1,0)/sqrt(D)=R_L/(2pi)^d.                       (6)

No factor1/2 is inserted. In particular opposite global signs cannot be
blindly paired: they generally select opposite sides of the fixed v.
The computation explicitly includes all64 signs.

## 4. Exact generic direction and orbit classification

Every integer fan ray has coordinate l1 norm at most12, checked against
the original rational rays. Thus each normal entry is at most12 in
absolute value. Two-column cross minors have absolute value at most288.
Use v=(1), (1009,1), or(1009^2,1009,1). A nonzero integer hyperplane normal
with coefficients bounded by288 cannot annihilate the last vector: its
leading nonzero base1009 coefficient dominates all lower terms. The
corresponding rank2 argument uses the smaller bound12. In addition to
this proof, every actual independent column facet in the enumeration is
formed exactly and its dot product with v is checked nonzero. Primitive
facet normals are retained in results.json.

Rank1 has one body-direction orbit. For rank2, flips of row representatives
replace Hamming k by6-k, leaving k=1,2,3. For rank3, normalize each physical
coordinate column by its first row sign; its projective class is one of
(1,1,1),(1,1,-1),(1,-1,1),(1,-1,-1). Signed permutations of the three
normal rows induce every permutation of these four classes (24 actions,
checked exactly). Physical coordinate permutations then identify triples
precisely by the sorted four class counts. Independence requires at least
three nonempty classes. The partitions of6 are exactly

    (4,1,1,0), (3,2,1,0), (2,2,2,0), (3,1,1,1), (2,2,1,1).

All4960 triples from the32 unoriented body representatives were separately
classified, and their exact Gram determinants checked positive. Their
orbit counts are480,1920,480,640,1440 respectively. The Gram determinants
are64,96,128,160,192. Thus no triple geometry was omitted.

For each representative, enumerate unique coordinate permutations of its
six projective columns, all64 signs, and21 cones. Every coordinate pattern
is weighted by product_class(count!) to restore all720 permutations. The
result is967680 logical terms for each of nine cases, with701568 unique
pattern/cone items in total. Rank1 multiplicity is720; rank2 multiplicity
is k!(6-k)!; rank3 uses the same class-factorial rule. Row permutations or
row-sign changes serve only the orbit classification and are not extra
multipliers in(5).

## 5. Certificate scope

The sole full enumeration checks the fan conversion, every body row's
constant sign on each cone, all triple classes, every normal facet and
all trace contributions. Only260 terms survive(5), and every survivor is
retained with its permutation/sign/cone identifiers, nonzero column set,
determinant, Cramer numerators and fully weighted exact rational value.
The minimal checker uses a separate SymPy matrix inverse/determinant and
fmpq summation to recheck these260 terms. It is not a second exhaustive
replay and does not by itself certify that no omitted term was positive;
that responsibility remains with the full enumeration source, count ledger
and this proof. Nine Arb512 comparisons with the stored B008 peak bounds
are postprocessing only. No physical Green/H kernel, tube grid or C6
estimate is evaluated or inferred.
