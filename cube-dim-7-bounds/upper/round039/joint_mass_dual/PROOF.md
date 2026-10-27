# R039 · Common-witness mass duality and geometric cluster payment

This is a general finite-family payment theorem and a precise route beyond separate pointwise load estimates. It does not claim a new uniform Cn bound. No production computation, n=7 tuning, MATLAB, or Git operation was performed. R037 remains frozen; its two clarifications were added separately as requested.

## 1. Fixed original history and its exact witness-mass problem

Let a>1 and let Qi=Q(yi,ri), 1≤i≤N, be selected small cubes from one actual history. Retain the actual disjoint increments

\[
I_i=Q_i\setminus\bigcup_{j<i}Q_j,\quad v_i=|I_i|,
\quad U=\sum_i v_i,\quad S_i=aQ_i,\quad V_i=|S_i|=(2ar_i)^n.
\]

Suppose the original nonnegative input f witnesses every Si at level λ>0, so ∫Si f≥λVi. Strict inequalities from the original maximal operator imply these weak inequalities; no new input class is assumed.

Define the exact common-witness mass cost

\[
\mathfrak m(\Gamma)=\min\left\{\int g:
g\in L^1_+,\ \int_{S_i}g\ge V_i\ \forall i\right\}. \tag{1}
\]

Then

\[
\frac{\lambda U}{\|f\|_1}\le\frac{U}{\mathfrak m(\Gamma)}. \tag{2}
\]

This is a selected-union payment ratio, not by itself a lower bound for Cn or a claim that the entire selected union belongs to the level set.

### Exact finite LP and dual

Partition a bounded box containing ∪Si by all coordinate face hyperplanes of the Si. Apart from a null set, this gives finitely many positive-volume open rectangular cells C. Every Si indicator is constant on each cell. Formula (1) is exactly the finite LP

\[
\min\sum_C m_C,\quad
m_C\ge0,\quad\sum_{C\subset S_i}m_C\ge V_i.
\]

The correspondence is exact: an L1 function gives its cell masses; conversely make g constant with value mC/|C| on each cell. There is a finite feasible solution, the objective is bounded below and its sublevel sets are compact in the nonnegative orthant, so an optimum exists. Finite LP strong duality gives

\[
\mathfrak m(\Gamma)=
\max\left\{\sum_i z_iV_i:z_i\ge0,
\ \sum_i z_i\mathbf1_{S_i}(x)\le1\ \text{a.e. }x\right\}. \tag{3}
\]

The same value is the infimum for strict witnesses, by multiplying any optimum by 1+ε. Thus (3) preserves one common f and all geometric overlap cells; it does not assign independent witness masses to cubes.

The old pointwise proof uses only the particular dual vector

\[
z_i=\frac{v_i/V_i}{\|W\|_\infty},\qquad
W(x)=\sum_i(v_i/V_i)\mathbf1_{S_i}(x).
\]

Its dual value is U/||W||∞. Optimizing the joint mass dual can be much stronger than this particular normalization.

## 2. A concrete geometric saving: an entire common-star cluster pays once

**Theorem 1 (uniform common-star payment).** If ∩iSi is nonempty, let R=maxiri and let j be any index of that radius. Then

\[
\mathfrak m(\Gamma)=V_j,\qquad
U\le\left(\frac{a+1}{a}\right)^nV_j,
\]

and every common witness f satisfies

\[
\lambda U\le\left(\frac{a+1}{a}\right)^n
\int_{S_j}f
\le\left(\frac{a+1}{a}\right)^n\|f\|_1. \tag{4}
\]

**Proof.** Pick x in the common open intersection. For every y∈Qi,

\[
\|y-x\|_\infty<r_i+ar_i\le(a+1)R.
\]

All actual increments therefore lie in Q(x,(a+1)R), and their disjointness gives U≤[2(a+1)R]^n. The j-th witness gives ∫Sj f≥λVj, proving (4). The mass LP costs at least Vj by that single constraint. Conversely concentrate an L1 density of mass Vj inside the nonempty common open intersection; it meets every constraint because Vi≤Vj. Thus the cost equals Vj. ∎

This theorem applies to the **entire** R037 high-load construction, since all stars share 0. It is stronger than its fixed-level spike diagnostic: every valid common witness pays the selected union with factor ((a+1)/a)^n, tending to 1 for fixed n as a→∞, even though W(0) can be (1−o(1))nlogn. Hence the loss is not forced by common-witness consistency; it is introduced by normalizing one pointwise load bound before using the volume geometry.

It still does not bound the whole maximal level set for an arbitrary f by 1: a general selected family need not have a single common-star intersection. The global shared mass must be accounted for when many clusters are used.

## 3. Fractional common-star clusters with paid representative overlap

Call J⊂{1,...,N} a cluster when ∩i∈JSi≠∅. Let j(J) be a chosen maximum-radius member, and put ca=((a+1)/a)^n.

**Theorem 2 (geometric cluster-payment certificate).** Suppose finitely many clusters J have coefficients θJ≥0 satisfying

\[
\sum_{J\ni i}\theta_J\ge1\quad\text{for every }i, \tag{5}
\]

and their chosen representative stars satisfy

\[
\sum_J\theta_J\mathbf1_{S_{j(J)}}(x)\le K
\quad\text{for a.e. }x. \tag{6}
\]

Then every original common witness f obeys

\[
\lambda\sum_i v_i\le c_aK\|f\|_1. \tag{7}
\]

**Proof.** The actual increments in any J remain disjoint, even though J is not an initial segment of the history. The geometric proof of Theorem 1 gives

\[
\sum_{i\in J}v_i\le c_aV_{j(J)}.
\]

Thus (5), the original witnesses and (6) imply

\[
\lambda U
\le\lambda\sum_J\theta_J\sum_{i\in J}v_i
\le c_a\sum_J\theta_J\lambda V_{j(J)}
\le c_a\int f(x)\sum_J\theta_J\mathbf1_{S_{j(J)}}(x)dx
\le c_aK\|f\|_1.
\]

No history is changed and no actual increment is replaced by an independently optimized set. ∎

The representative-overlap factor K is the explicit cost of reusing one mass source. This is a joint geometric payment certificate, rather than another pointwise bound for the old W. A certificate with K=o(nlogn) uniformly over the relevant selected histories would improve the framework; such a uniform construction is not proved here.

For the original threshold covering algorithm, the already accepted cover step gives

\[
C_n\le\left(1+\frac{(\beta+1/a)^n}{\tau}\right)c_a K_*, \tag{8}
\]

provided K* is a uniform valid certificate bound over all its selected histories. The unknown supremum K* must not be replaced by values from a finite instance search.

## 4. A finite, geometry-exact LP interface without free cluster addition

The axis-parallel star family has Helly number 2: a finite subfamily intersects iff each pair intersects. For open cubes the pair test is ||yi−yj||∞<a(ri+rj). The strict pairwise interval intersections give a nonempty open intersection in each coordinate.

For a complete finite certificate interface, use the same positive-volume arrangement cells C from Section 1. Let TC={i:C⊂Si}. For every cell C and every j∈TC define the radius-truncated cluster

\[
J_{C,j}=\{i\in T_C:r_i\le r_j\}.
\]

It has a common pivot in C and chosen representative j. Any admissible cluster J with representative j can be enlarged to one of these clusters **without changing j**: choose a pivot in its common open intersection away from all finitely many face hyperplanes. The resulting J_C,j contains J. Thus restricting to these finite columns loses no certificate of the form in Theorem 2; restricting just to arbitrary maximal cliques without preserving the representative would not have this property.

The exact cluster LP is

\[
\min K,\quad\theta_{C,j}\ge0,
\]

\[
\sum_{(C,j):i\in J_{C,j}}\theta_{C,j}\ge1
\quad\forall i,
\]

\[
\sum_{(C,j):j\in T_D}\theta_{C,j}\le K
\quad\forall\text{ arrangement cells }D.
\]

The last constraints explicitly charge representative overlap. One feasible nonnegative certificate is enough for an upper bound; no claim of a global optimizer is needed. Floating-point output would need outward checking of these linear inequalities. This report provides no production solve.

An equivalent continuous notation uses pivot measures θj supported in Sj: item i is covered by Σj:rj≥ri θj(Si∩Sj)≥1, and the paid load is Σj θj(Sj)1Sj≤K. The finite cell form above avoids measurability or selector ambiguity.

## 5. Counterexample diagnostic for unpaid representative mass

Take two disjoint equal small cubes with centers ±2r e1 and radius r, and choose a>2. Their two star cubes overlap in a nonempty open set. An L1 spike in that intersection, of mass (1+ε)V with V=(2ar)^n and ε>0, strictly witnesses both at level 1. Nevertheless

\[
\int_{S_1}f+\int_{S_2}f=2\|f\|_1.
\]

Hence two separately paid representative constraints cannot be added and then charged once to ||f||1 without an overlap factor. This example diagnoses the invalid mass-charging step; it is not a Cn lower-bound construction. Treating the two items as a single common-star cluster gives the correct K=1 certificate and illustrates why cluster compression and global overlap payment must be performed together.

### A center-pivot greedy cluster rule can lose an arbitrary factor

The pivot need not be one of the original centers. This matters even for legal selected histories in every dimension. Fix N, set Ri=3^(-i), yi=(3/4)Ri e1, ri=Ri/a, and choose a>max(3,2t) with (N−1)a^-n≤τ. The small Qi are disjoint: in their first coordinate,

\[
(3/4+1/a)R_{i+1}<(3/4-1/a)R_i
\]

whenever a>8/3. Their actual weights are a^-n. Every star Si has radius Ri and contains 0. Each earlier bQi also contains every later center, because ||yi−yj||<(3/4)Ri<βRi. Therefore every item remains eligible and the actual largest-radius algorithm selects all N.

Consider instead the proposed cluster rule “take the largest remaining star, use its own center as pivot, and remove all stars containing that pivot.” No later Sj contains yi: its first-coordinate right endpoint is at most (7/4)Rj≤(7/12)Ri<(3/4)Ri. The rule therefore creates N singleton clusters, whose representative overlap at 0 is N. Yet the true common pivot 0 gives a single cluster and K=1.

This is a genuine geometric/history counterexample to that center-only clustering strategy, not merely an abstract overlap diagram. It does not contradict Theorem 2: the poor strategy has paid K=N. Radius order and eligibility alone do not make this particular pivot strategy efficient. The radius-truncated arrangement-cell interface retains the missing common pivot.

## 6. Scope and next unresolved step

Theorems 1–2 are valid for every n, every finite actual history, all positive radii, arbitrary L1 witnesses, and every a>1. They do not use a bounded-radii assumption; a common homothety preserves all statements. Eligibility has not yet been exploited in the cluster LP, so the theorem is a sound general interface rather than a claimed improvement for all selected histories.

The next substantive task is to construct geometry- and history-dependent pivot measures or cluster coefficients whose representative overlap is smaller than Θ(nlogn), while preserving every item and accounting for the common f exactly. This could use the actual radius order and rejection information. Merely minimizing the old W, or adding the common-star payments without (6), cannot provide that conclusion.

