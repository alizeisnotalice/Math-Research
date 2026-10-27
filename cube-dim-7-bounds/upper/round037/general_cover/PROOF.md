# R037 · Real-history high-shell saturation and a pointwise-load scope barrier

Proposed self-contained result for independent review. This is a barrier for the current selection algorithm's **universal pointwise weighted load estimate**, not a lower bound for the maximal-operator constant Cn. It retains actual increments, nonincreasing radii, all intermediate eligibility tests, a common load support point, and a genuine positive-weight final minimum-radius cube. No n=7 work, production computation, MATLAB, or Git operation.

## 1. The exact statement

Fix integer n≥3, t>S≥2, β=1+1/t, and threshold τ>0. Put

\[
m(z)=\max(1,\|z\|_\infty),\qquad
D=\{\|z\|_\infty<t:\|\mathbf1-z\|_\infty>\beta m(z)\},
\]

\[
d\nu(z)=2^{-n}m(z)^{-n}
\mathbf1_{D\cap\{S<m(z)<t\}}\,dz,
\quad M=\nu(\mathbb R^n),\quad
M_0=n\log(t/S),
\]

\[
\Gamma=\frac{(n-2)S}{16t},\qquad
\varepsilon_{\rm bad}=\frac{4M_0^2e^{-\Gamma}}{\tau}.
\]

Use exactly the original finite selection algorithm: a>2t, b=aβ, actual increments Ij=Qj\(∪i<jQi), weights wj=|Ij|/(2arj)^n, and eligibility Σi<j wi 1_bQi(yj)≤τ, with largest eligible radius selected first.

**Theorem 1 (actual finite histories).** For each ε>0 there is a0 such that every real a≥a0 admits a finite original input family Q1,...,QN,Qk for which:

1. The algorithm selects the entire family in nonincreasing radius order, with rk=1 the unique minimum and rj>S for j≤N.
2. All small cubes are pairwise disjoint, so their actual increments satisfy Ij=Qj and every weight is exactly a^-n.
3. The point 0 belongs to every aQj, including aQk. The final item has genuine positive weight a^-n.
4. The final center yk is uncaptured by every earlier bQj, so its eligibility sum is zero.
5. The actual pointwise weighted load obeys

\[
W(0)=\sum_jw_j\mathbf1_{aQ_j}(0)
\ge M-\varepsilon_{\rm bad}-\varepsilon
\ge B_n(t)-1-n\log S-\varepsilon_{\rm bad}-\varepsilon. \tag{1}
\]

Thus the continuum model below is restored to real finite histories, not left as an abstract packing relaxation.

## 2. A full incoming-energy upper bound

Define the full incoming potential of the base cap measure by

\[
T(w)=\int\mathbf1_{\{\|z-w\|_\infty<\beta m(z)\}}\,d\nu(z).
\]

No radial order is imposed in T. Bounding the entire incoming load will automatically bound every previous-only eligibility sum.

**Lemma 2.**

\[
E:=\int T(w)d\nu(w)\le M_0^2e^{-\Gamma}. \tag{2}
\]

**Proof.** Enlarge both measures to the full-shell measure 2^-n||z||^-n dz on S<||z||<t, whose mass is M0. Conditional on radii s≤ρ and on the two chosen maximal faces, all but at most two coordinate pairs are independent uniforms Z∼Uniform[−ρ,ρ], W∼Uniform[−s,s]. If the directed source has the smaller radius, its true capture event is contained in the enlarged event with radius βρ, so it suffices to use β times the larger radius in either orientation.

For each remaining free coordinate pair, two exact tail triangles give

\[
\mathbb P(|Z-W|\ge\beta\rho)
=\frac{(s-\rho/t)^2}{4s\rho}.
\]

Here ρ/t≤1≤s/2 because s≥S≥2. Thus the failure probability is at least s/(16ρ), and the success probability is at most exp[−s/(16ρ)]. Full capture implies success in the n−2 free pairs. Its conditional probability is at most

\[
\exp[-(n-2)s/(16\rho)]\le e^{-\Gamma}.
\]

Integrating both full radial masses proves (2). Restriction to the cap only decreases the integral. The prefactor is M0², not M²: the cap restriction does not preserve independent coordinate laws. ∎

Markov now gives

\[
\nu\{T\ge\tau/4\}\le4E/\tau\le\varepsilon_{\rm bad}. \tag{3}
\]

Moreover T is continuous on Rn. On each max-coordinate cone of z, every capture boundary is contained in one of finitely many nonzero-normal hyperplanes

\[
z_i-w_i=\pm\beta(\pm z_j).
\]

When i=j, the coefficient has magnitude at least β−1>0; when i≠j, it is also nonzero. Each such boundary is Lebesgue-null. Dominated convergence proves continuity. Therefore

\[
G=D\cap\{S<\|w\|_\infty<t\}\cap\{T(w)<\tau/4\}
\]

is open and ν(G)≥M−εbad.

## 3. Compact blocks and true disjoint increments

The following construction proves Theorem 1. All choices before a are finite and made at fixed n,t,S,τ,ε.

Choose a finite union K of closed axis-parallel blocks with disjoint interiors and small gaps, compactly contained in G, such that ν(K)≥ν(G)−ε/4. This is possible by inner regularity and box approximation of the open set G. If the desired right side of (1) is nonpositive, a single final cube already proves the assertion; below consider the nontrivial case.

Refine the blocks sufficiently. To block Pl assign a constant radius Rl satisfying

\[
R_l>\sup_{p\in P_l}\|p\|_\infty,
\quad S<R_l<t,
\quad 0<R_l-m(p)\le\delta\quad(p\in P_l),
\]

where δ>0 can be chosen arbitrarily small after refinement. Because K is compactly inside the strict cap D, choose a fixed η>0 and uη=(1−η)1 such that, for all blocks and p∈Pl,

\[
\|u_\eta-p\|_\infty>\beta R_l. \tag{4}
\]

Indeed the minimum of ||1−p||−βm(p) on K is positive; take η+βδ below that margin.

The continuous model associated with the blocks is

\[
T_R(w)=\sum_l\int_{P_l}
\frac{\mathbf1_{\{\|p-w\|_\infty<\beta R_l\}}}{(2R_l)^n}\,dp.
\]

Since Rl≥m(p), its density is no larger than the base density. Inflating capture radii from βm(p) to βRl only adds a uniformly small boundary layer. More explicitly, for all w the extra region is contained in a union of at most 4n² slabs of width O(δ/(β−1)) about the nonzero-normal hyperplanes from Section 2, inside [−t,t]^n. Its weighted measure is ≤C(n,t,S,β)δ. Thus

\[
T_R(w)\le T(w)+C\delta\quad\text{uniformly in }w. \tag{5}
\]

Choose δ small enough that Cδ<τ/16 and that the block mass Σ|Pl|/(2Rl)^n≥ν(K)−ε/4. On K, (5) gives TR<5τ/16<τ/2.

For each a, tile each Pl with p-space cubes of side 2Rl/a and keep only cells fully contained in Pl. If pi is a retained center, define the original small cube

\[
Q_i=Q(ap_i,R_l).
\]

Different retained cells give disjoint original open Qi. Consequently their true increments are Qi themselves and wi=a^-n. Also ||pi||<Rl implies 0∈aQi, and (4) implies a uη∉bQi.

### Uniform discrete eligibility, including all source items

For any center w and every retained source cell Ci in block Pl, if ||w−pi||<βRl, then every p∈Ci satisfies ||w−p||<βRl+Rl/a. Therefore

\[
\sum_{i:\|w-p_i\|<\beta R_i}a^{-n}
\le\sum_l\int_{P_l}
\frac{\mathbf1_{\{\|w-p\|<\beta R_l+R_l/a\}}}{(2R_l)^n}\,dp.
\]

The difference between this upper bound and TR(w) is O(1/a), uniformly in w: each block contributes only a boundary strip of a fixed axis-parallel box. For all sufficiently large a, the entire source load at every retained center w∈K is below τ/2, hence below τ. This includes all sources, even itself and later sources, and is stronger than the algorithm's actual previous-only condition.

Thus every remaining source cube stays eligible throughout selection. The largest-radius rule selects them in nonincreasing radii. Append the final cube

\[
Q_k=Q(a u_\eta,1).
\]

For large a it is disjoint from all source cubes, since K is a positive distance from uη. It has 0∈aQk, positive weight a^-n, and unique minimum radius. Every source radius is >S≥2, so all sources are selected first. At its selection the final cube has zero captured sum by (4), and hence is selected too.

Finally the retained grid-cell mass tends to Σ|Pl|/(2Rl)^n as a→∞, with only boundary loss. For all sufficiently large a it is at least ν(K)−ε/2≥M−εbad−ε. This proves the first lower bound of (1). The omitted unit core has mass≤1 and the omitted shells 1<m≤S have mass≤nlogS; therefore M≥Bn(t)−1−nlogS, proving the second bound. ∎

The quantifier order is important: first fix n,t,S,τ,ε, then choose the compact blocks, radii, η and boundary margins; finally choose a0. For every a≥a0 a finite family is obtained. No optimizer or selector convergence is assumed.

## 4. A uniform near-linear-buffer scope barrier

Let L=logn and fix C≥0. Suppose

\[
\frac n{L^C}\le t\le nL^C,\qquad\tau\ge1.
\]

Take S=L^(C+3). For all sufficiently large n, 2<S<t and

\[
\Gamma\ge cL^3,\qquad M_0\le2nL,
\qquad\varepsilon_{\rm bad}=o(1)
\]

uniformly in this whole t window and all τ≥1. The accepted integral bound

\[
I_n(t)\le\log_+(t/n)+8+\log2
\]

implies Bn(t)=nL+O_C(nlogL) throughout this window. Combining with (1), real finite histories satisfy

\[
W(0)\ge n\log n-O_C(n\log\log n)-o(1)-\varepsilon. \tag{6}
\]

In particular at the standard t=n/(L−3logL), τ=√(AnBn), choosing S=L² already suffices: Γ≥cL³ and

\[
W(0)\ge B_n(t)-2n\log L-1-o(1)-\varepsilon.
\]

The original universal pointwise upper bound tends to τ+Bn(t), which is (1+o(1))nlogn at that standard parameter. Consequently its leading pointwise load is genuinely saturated by actual histories.

**Precisely what is excluded.** A universal pointwise W estimate for this same algorithm, claimed for every original history and every positive common scale, cannot reduce a fixed positive proportion of nlogn in the displayed near-linear t window. Actual increment disjointness, center-radius coupling, the final support item, and all intermediate eligibility conditions are already present in the construction.

This does **not** exclude improved simultaneous use of the covering measure and f-witness inequalities, a different selection algorithm or weighting, or t regimes outside the stated window. It is not a lower bound for Cn. The original maximal operator still ranges over every positive radius; constructing finite legal histories in a particular radius band is enough to test a purported universal pointwise bound. A common homothety gives any desired final positive radius.

## 5. Counterexample diagnostic: large W is compatible with a near-unit weak ratio

Every constructed aQi contains 0 in its interior. A sufficiently small L1 spike f supported near 0 can simultaneously witness all these finite star cubes: choose its mass M_f>λ(2aRmax)^n, where Rmax is the largest selected radius. This is a genuine common f, not independent mass assignments to individual cubes.

Nevertheless a small-support spike has weak ratio at most 1+o(1). If f≥0 has mass Mf and support in (−εf,εf)^n, then M_n f(x)>λ implies

\[
\|x\|_\infty<\tfrac12(M_f/\lambda)^{1/n}+\varepsilon_f.
\]

Indeed any witnessing cube must have radius at least ||x||−εf while its average is at most Mf/(2R)^n. Hence

\[
\frac{\lambda|\{M_nf>\lambda\}|}{M_f}
\le\left[1+2\varepsilon_f(\lambda/M_f)^{1/n}\right]^n
\longrightarrow1\quad(\varepsilon_f\downarrow0).
\]

Thus the very same high-load history admits legitimate common witnesses whose actual weak ratio is near 1. This does not prove every witness for that history is harmless; it diagnoses the invalid inference “large pointwise W gives a comparable Cn lower bound.” The finite proof pays all scales at the single location 0, whereas the f-witness/covering efficiency can have a different scale balance. That joint balance, or a genuinely different cover selection, is the remaining route beyond this scope barrier.

## 6. Status and reproduction interface

All statements above are analytic; no new production computation was run. The only needed computational interface, if the main agent requests diagnostics, is the exact one-coordinate identity in Section 2 and the explicit quantities M0, Γ, εbad. No old grid or low-dimensional tuning is needed. All new files are restricted to `upper/round037/general_cover/`; R030–R036 remain unchanged.

