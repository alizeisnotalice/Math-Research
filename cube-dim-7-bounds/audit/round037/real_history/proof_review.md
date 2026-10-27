# R037 real-history load barrier: independent full-scope audit

2026-10-28. **ACCEPT as an actual finite-history pointwise-load barrier.**
The continuum model is restored to legal finite instances of the original
selection algorithm; the conclusion need not be weakened to a relaxation
barrier. The theorem does **not** give a lower bound for C_n, and does not
exclude joint covering/witness improvements, different algorithms or
weights, other parameter regimes, or lower-order pointwise improvements.

## 1. Source and precise accepted quantifiers

Source read in full:
`/Users/zhengzhihao/Desktop/T/cube-r030-geometry/cube-dim-7-bounds/upper/round037/general_cover/PROOF.md`,
SHA256 `6261341afbdafe868616e3e365f3484e0dab7e9e8e6ca09f83c56c93751df90e`.

For each fixed integer n>=3, t>S>=2, tau>0 and epsilon>0, there is a finite
a0 such that **every real a>=a0** admits the stated finite history. All
source radii exceed S; the last radius 1 is uniquely minimal and has
positive weight. All increments are the actual disjoint cubes; all star
cubes contain 0; every intermediate selection is eligible; and the final
center is uncaptured by all preceding source cubes.

For fixed C>=0 and the near-linear window n/L^C<=t<=nL^C, tau>=1, the
large-n constants in the leading load bound and the bad-mass estimate are
uniform in that window. The finite construction may have
a0=a0(n,t,S,tau,epsilon); a common a0 for an entire parameter window is
neither needed nor asserted. No upper bound on the required size of a0
or number of cubes is proved.

Only `audit/round037/real_history/` is written. The audit is analytical;
no production calculation, n=7 work, parameter search, MATLAB or Git
operation was performed. Frozen R035/R036 materials were not changed.

## 2. Full incoming-energy estimate is valid

The full weighted shell measure on S<||z||_infty<t has mass
M0=n log(t/S). Conditional on the two radii s<=rho and the two maximal
faces, all but at most two coordinate pairs remain independent uniform
variables on [-rho,rho] and [-s,s]. The conditioning is legitimate for
this enlarged full-shell law, not for the cap-restricted law.

If the directed source has the smaller radius, replacing its capture
radius by beta*rho only enlarges its event. In a free coordinate, put
v=rho/t. Since v<=1<=s/2, the two tail triangles have exact failure
probability

```
Pr(|Z-W|>=rho+v)=(s-v)^2/(4s*rho)>=s/(16rho).
```

The triangles fit inside the rectangle because rho>=s and 0<=v<s.
The two areas sum before division by its total area 4s*rho, so there
is no missing factor 2. Capture requires success in the n-2 free
pairs, and hence has conditional probability at most
exp[-(n-2)s/(16rho)]<=exp(-Gamma). Integrating over *all ordered pairs*
of full shells gives M0^2 exp(-Gamma). No additional ordering factor is
needed; the radius ordering was only a conditional argument. Restriction
back to the cap only decreases this nonnegative integral.

Markov therefore removes at most epsilon_bad=4M0^2 exp(-Gamma)/tau
mass by requiring T(w)<tau/4. Its measure is correctly bounded with M0,
not with the smaller cap mass M; cap conditioning would not preserve
the independent coordinate law.

## 3. Openness, block labels and the uniform inflation bound

For each fixed w, the capture-boundary set in z is contained, on the
finitely many maximal-coordinate cones, in affine hyperplanes
z_i-w_i=+/- beta(+/- z_j). Each normal is nonzero. Even when i=j,
its coefficient magnitude is at least beta-1>0. The source measure is
absolutely continuous and finite, so dominated convergence proves T is
continuous. The retained G is consequently an open subset of the strict
cap and open radial band.

Finite closed axis-parallel blocks compactly inside G approximate its
weighted measure arbitrarily well. Small separating gaps can be arranged
with arbitrarily small loss. After a fixed compact union K is chosen,
refining it makes the oscillation of m on each block arbitrarily small.
One can therefore choose constant labels R_l with

```
sup_{p in P_l}m(p)<R_l<t,
0<R_l-m(p)<=delta on every block,
R_l>S.
```

Here delta may be chosen as small as needed. The strict cap margin on
compact K is positive. Taking eta+beta*delta below it yields
||u_eta-p||>beta R_l for all p in every block, with
u_eta=(1-eta)e inside the unit cube. This simultaneously preserves a
real final support point and its strict uncaptured margin.

The critical comparison T_R<=T+C delta is indeed uniform in w, with
C depending only on fixed n,t,S,beta rather than the number of refined
blocks. The block density (2R_l)^(-n) is no larger than the base density
(2m(p))^(-n). Extra capture from radius inflation is possible only when

```
0<=||p-w||_infty-beta*m(p)<beta*delta.
```

On maximal-coordinate and distance-coordinate cones this set lies in
at most 4n^2 slabs with width O(delta/(beta-1)), inside the fixed cube
[-t,t]^n. Their weighted volume is bounded uniformly in their affine
positions w_i. Using this single global slab union avoids a factor
depending on the number of blocks; there is no circular choice of delta
and block refinement. The constant can be very large as beta approaches
one, but n,t are fixed before delta and a are selected.

The block mass converges to the base K mass as delta decreases, by the
pointwise uniform proximity of R_l to m(p) on S<m<t. Thus both the load
error C delta<tau/16 and the stated mass error can be paid with one
sufficiently small delta.

## 4. Discrete restoration controls every source's eligibility

For each real a, tile each block using cells of side 2R_l/a and retain
those fully contained in it. The physical source cube is exactly a
times its p-space cell. Consequently all source cubes are disjoint open
sets even across blocks with different radius labels, and their actual
increments equal themselves under every ordering. Their weights are
exactly a^(-n).

At a normalized test center w, each captured retained source center p_i
satisfies that its *whole cell* lies in the box of radius
beta R_i+R_i/a about w. Its density integral over that cell equals

```
(2R_i)^(-n)*(2R_i/a)^n=a^(-n).
```

Summing cells therefore gives the source's upper integral with inflated
box radii. The difference from T_R is O(1/a) uniformly in w, because
there are now a fixed finite number of blocks and each contributes only
a boundary strip of a fixed box. All refinement choices preceded a.

Every retained center lies in K, where T<tau/4 and T_R<5tau/16. For all
sufficiently large a, the total discrete incoming source load is below
tau/2 at every retained center. Crucially, this includes itself and all
future sources. Every actual previous-only eligibility sum is a smaller
nonnegative subsum. Thus no source can become ineligible, and the
largest-eligible-radius algorithm selects every source in nonincreasing
label order, including arbitrary consistent tie choices. This is the
step that makes the construction an original history rather than a
packing relaxation.

The strict label inequality ||p_i||<R_i ensures 0 belongs to every aQ_i.
For every sufficiently large real a, the retained mass tends to
sum_l |P_l|/(2R_l)^n. This is ordinary boundary-cell loss for a mesh tending
to zero; it does not require a to follow an integer subsequence.

## 5. Final item, exact increments and error accounting

Append Q_k=Q(a u_eta,1). All source labels exceed S>=2, so its radius is
uniquely minimal and it is selected last once the sources are selected.
It has 0 in aQ_k because ||u_eta||<1. The fixed cap margin places its
center outside every bQ_i. Its eligibility sum is therefore zero at
every prefix, including at its own selection.

The fixed compact source union stays a positive distance from u_eta;
after physical scaling by a, the final cube is disjoint from every source
for all sufficiently large a. Hence its actual increment is also the
whole cube and its positive weight is a^(-n). No zero-weight or fictitious
terminal item is used.

The losses in passing from G to K, from base density to radius labels,
and from blocks to retained cells can each be made smaller than the
assigned epsilon fractions. Their sum is at most epsilon; adding the
positive final weight only helps. This proves

```
W(0)>=M-epsilon_bad-epsilon
    >=B_n(t)-1-n log S-epsilon_bad-epsilon.
```

The last inequality merely bounds the omitted unit core by 1 and the
shells between 1 and S by n log S. If the desired right side is
nonpositive, the separately stated single-final-cube case is legitimate.
In the large-n window the right side is positive and the source family
is nontrivial.

## 6. Near-linear-window quantifiers and the pointwise barrier

For fixed C>=0, choose S=L^(C+3). Eventually S<t uniformly in
n/L^C<=t<=nL^C, and

```
Gamma>=(1-2/n)L^3/16,
M0<=2nL,
epsilon_bad<=16n^2L^2 exp(-c L^3)=o(1)
```

uniformly for tau>=1. The large-n threshold may depend on C. The accepted
integral bound in the source gives B_n(t)=nL+O_C(n log L) throughout this
window: n log t differs from nL by at most Cn log L, and the subtracted
integral contributes at most another O_C(n log L). Also n log S is
(C+3)n log L. Thus the claimed real-history lower bound has the correct
uniform leading coefficient.

At the standard parameter t=n/(L-3log L), using S=L^2 yields the stronger
displayed loss 2n log L and Gamma of order L^3. The accepted limiting
pointwise budget tau+B_n(t) has leading term nL there, so this leading
pointwise load is saturated by actual histories. This does not contradict
the R033-R036 lower-order deficit improvements; the construction itself
still loses O(n log L), much larger than their certified corrections.

The obstruction concerns universal *pointwise W estimates for this same
algorithm* in the stated t window and tau>=1. It is not a theorem about
all scale-sensitive methods, thresholds below one, arbitrary t regimes,
all joint covering certificates, or the optimal operator constant.

## 7. A genuine common input and the spike diagnostic

There are finitely many open star cubes, and each contains 0 in its
interior. Their intersection therefore contains a small neighborhood of
0. A nonnegative L1 spike supported inside it, with mass greater than
lambda times the largest star volume, is simultaneously a strict witness
for every star cube. No independent or incompatible per-cube mass
assignments occur.

Put R0=(M_f/lambda)^(1/n)/2. If the spike support lies in
(-epsilon_f,epsilon_f)^n, any strict maximal witness has radius below
R0 and must meet that support. The maximal level set is consequently
contained in ||x||<R0+epsilon_f, proving the source's ratio upper bound
and its limit at most one. A matching lower inclusion
||x||<R0-epsilon_f follows by choosing a radius between
||x||+epsilon_f and R0 that contains all spike mass. Thus the ratio can
indeed be made to tend to one for the same constructed history.

For this diagnostic n and the finite history are fixed while the spike
support shrinks; no uniform-in-n spike size is asserted or needed. The
near-unit upper bound alone already prevents inferring a comparable C_n
lower bound from W(0). The source correctly does not claim that every
possible witness input for this history is harmless.

## 8. Final status

The continuum-to-finite-a grid, radius labels, all-source eligibility,
true final increment and common-input restoration are closed. The stated
result is an actual-history barrier with the scope above, not merely a
relaxation barrier. No numerical evidence or unproved convergence of
selectors was used. No new C_n upper bound or C_n lower bound follows
from this theorem by itself.
