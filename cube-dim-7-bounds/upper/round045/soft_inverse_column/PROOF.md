# R045 draft · Soft inverse-column weights with a recomputed row budget

Status: complete proposed full-operator mechanism, submitted for independent review. No numerical threshold is claimed. This proof does not replace the current frozen upper bound until reviewed. It uses a bounded position-dependent source weight, then recomputes both row and column scales. It does not retain a fixed fraction of the old unweighted energy. Therefore the R043 product-moment certificate barrier does not apply to this payment scheme.

No production computation, n=7 work, MATLAB, or Git operation. Other agents' independent sources are not edited.

## 1. Proposed conclusion and exact proof parameters

First fix an absolute constant

\[
0<d\le2^{-24},\qquad\varepsilon_0=2^{-20}.
\]

As n→∞ put

\[
L=\log n,\quad x=(1-32d^2)L,\quad t=n/x,
\quad\beta=1+1/t,\quad\zeta=n^{-3},\quad\kappa=\beta-\zeta,
\]

\[
h=\frac{\varepsilon_0}{Ld},\qquad
\theta=e^{-12Ld^2},\qquad k=nh\theta,
\quad I=[d/4,3d/4].
\]

Let A0=β^n and B0=B_n(t) be the accepted original cap budgets. The proposed result is, for every sufficiently large n depending only on the fixed d,

\[
C_n\le F\left(A_0,B_0-\frac{nd}{65536}\right),
\qquad F(A,B)=(\sqrt A+\sqrt B)^2. \tag{1}
\]

In particular choosing the explicit absolute d=2^-30 would give, for all sufficiently large n,

\[
C_n\le n\log n-n\log\log n-2^{-47}n. \tag{2}
\]

This is a candidate until independently accepted. The common sufficiently-large threshold can be enormous and has not been quantified. No claim is made that (2) dominates all previously valid bounds at their original finite thresholds.

## 2. General weighted energy and the unchanged actual row payment

Let μ≤λ be the actual-increment and ambient measures on one common history, and let A,P be disjoint. Suppose a raw Borel kernel K0∈[0,1] is supported in A×P and has actual row bound ∫K0(w,z)dμ(z)≤τ for μ-a.e. w∈A.

For Borel row and source weights 0≤a(w),b(z)≤1 use

\[
K(w,z)=a(w)b(z)K_0(w,z).
\]

Then actual row payment is at most a(w)τ; it is not multiplied by an unbounded inverse-column factor. If Rλ,Cλ are degrees of K and both are at most M, the accepted deficit proof gives

\[
\Delta\ge\frac{E_{\lambda\lambda}-\tau\int_Aa\,d\lambda}{M}. \tag{3}
\]

Indeed Eμμ≤τ∫a dμ, and Eλλ−Eμμ≤M[d(A)+d(P)]≤MΔ for the same missing measure d=λ−μ. More precisely the paid deficit potential is a(w)(R_b(w)−τ) on A and b(z)C_a(z) on P. The application below uses a=1 and a bounded source weight b.

## 3. Ambient geometry and the soft source weight

Use the accepted common corner measure

\[
d\lambda(z)=2^{-n}m(z)^{-n}\mathbf1_{\mathcal D}(z)dz,
\quad m(z)=\max(1,\|z\|_\infty),
\]

\[
\mathcal D=\{\|z\|_\infty<t:\exists i,z_i<1-\beta m(z)\}.
\]

Before cap restriction each face at radius R>1 has measure dR/(2R), with n−1 independent Uniform[−R,R] free coordinates.

Targets have r=||w||∞=1+u, u∈I. Sources have ρ=||z||∞=tq, q=1−v, v∈I. Radial matching is

\[
|u+v-d|<h,\quad\text{equivalently }|r-q-d|<h. \tag{4}
\]

For every source radius set

\[
R(q)=q+d-h,\quad
\mathfrak a(q)=(\kappa-1)\rho=q(1-\zeta t),\quad
g(q)=R-\mathfrak a=d-h+\zeta tq.
\]

Every matched target radius is at least R and at most R+2h. Eventually R>1, ρ≥4r, 0<g≤d, h+ζt<d/2, q>1/2 and κ>1. Define exact normalized target-coordinate column fractions

\[
\ell_R(z_i)=\min\left(1,\frac{R+\kappa\rho-|z_i|}{2R}\right),
\quad G(z)=\prod_i\ell_R(z_i)>0,
\]

and the bounded soft weight

\[
b(z)=\min(1,\theta/G(z)). \tag{5}
\]

The full target-column estimate from the interval width and deleted-coordinate products is

\[
C_{\rm raw}(z)\le\frac{4nh}{R}G(z).
\]

Here normalized lengths decrease as target r increases, each length at R is ≥1/2, and the matched target log-width is ≤2h/R. Hence (5) yields for every source point

\[
b(z)C_{\rm raw}(z)\le4nh\theta=4k. \tag{6}
\]

This bounds columns without a uniform high-column truncation H and without assigning an independent mass source to each target.

## 4. Exact one-coordinate inverse row

For fixed q,R,ρ,g define

\[
J(w;q)=\frac1{2\rho}\int_{[-\rho,\rho]\cap(w-\kappa\rho,w+\kappa\rho)}
\frac{dz}{\ell_R(z)}.
\]

Write

\[
A(R,g,T)=
\begin{cases}
2R\log(1+T/(2R-g)),&0\le T\le g,\\
2R\log(2R/(2R-g))+T-g,&T\ge g,
\end{cases}
\]

\[
B(R,g)=A(R,g,g)-g,\qquad
T=(|w|-\mathfrak a)_+.
\]

At most one source edge is removed because ρ≫r and κ>1. The exact formula is

\[
\boxed{J(w;q)=1+\frac{2B(R,g)-A(R,g,T)}{2\rho}.} \tag{7}
\]

It is equivalent to the stable closed form with the two one-sided excesses in the independent R046 diagnostics. Formula (7), not a sampled numerical fit, is used below.

### 4.1 Baseline at target radius r=R

For W uniform on [−R,R], Fubini gives exactly EJ(W;q)=1: for a fixed z, its fraction of accepted W is precisely ℓ_R(z), canceling the inverse weight.

When g≤R/2 and ρ≥4R, B≤g²/(3R), A(R,g,T)≤(4/3)T for 0≤T≤g, and |J−1|≤g/ρ≤1/8. Directly integrating the tail T=(|W|−(R−g))_+ gives a conservative bound

\[
\mathbb E(J-1)^2\le3g^3/\rho^2\quad(R\ge1).
\]

For example use (B−A/2)^2≤2B²+A²/2 and ET²=g³/(3R). Thus the elementary log bounds imply

\[
-6g^3/\rho^2\le\mathbb E\log J\le0,
\qquad\operatorname{Var}(\log J)\le12g^3/\rho^2. \tag{8}
\]

At the pinned coordinate W=±R, logJ=O(g/ρ). These conditional-face bounds retain the fixed coordinate and the nonzero mean of logJ.

### 4.2 Uniform matching-radius bridge, with kinks handled by absolute continuity

Fix target r and |w|≤r, and let q range between r−d−h and q0=r−d+h. At q0, R=r. For the entire interval, R∈[r−2h,r], R,q≥1/2, g≤d, ρ≥4R, and T≤2d. Except at the capture-onset kink |w|=mathfrak a(q), differentiating (7) gives the conservative bound

\[
\left|\partial_q\log J(w;q)\right|
\le\frac{64}{\rho}
\left[d^2+\mathbf1_{\{|w|>r-2d\}}\right]. \tag{9}
\]

The function is locally Lipschitz, so (9) holds a.e. and may be integrated across its kinks. Here is one explicit coarse derivation. For f(y)=2R/(2R−g+y), 0≤y≤g, use |f_R|≤8d, |f_g|≤4, f≤2. The piecewise A has |A_R|≤16d², |A_g|≤8d, A_T≤2; B has |B_R|≤8d², |B_g|≤4d. Since R'=1, g'=ζt≤d and |T'|≤1{|w|>r−2d}, one gets |J'|≤32[d²+1{tail}]/ρ. Also J≥1/2. This proves (9). At T=g the two branches have matching derivatives; the onset T=0 is handled a.e.

## 5. Typical targets have uniformly bounded full inverse rows

For a target w with r=||w||∞ put q0=r−d+h and evaluate each J(wi;q0) with R=r. Define

\[
\mathcal A=\mathcal D\cap\{u\in I,
|\sum_i\log J(w_i;q_0)|\le1/8,
\quad K_w:=\#\{|w_i|>r-2d\}\le4nd\}.
\]

For each maximal target face, the n−1 free coordinates are independent Uniform[−r,r], while the pinned one equals ±r. By (8), both the absolute mean and variance of their logarithmic sum are O(nd³/t²)=O_d(L²/n)→0, uniformly over the compact radial interval. The pinned term tends to zero. Chebyshev therefore gives the baseline logarithmic condition with probability 1−o(1).

Moreover Kw=1+Bin(n−1,2d/r), with mean≤1+2nd and variance≤2nd. At threshold 4nd, its bad probability tends to zero. The target corner-captured part is exponentially small in n. Thus eventually

\[
a_*:=\lambda(\mathcal A)\ge nd/16. \tag{10}
\]

Every matched q lies between q0−2h and q0. Sum (9) over coordinates, use Kw≤4nd and n/ρ≤2L, and integrate over a q-length≤2h. This gives

\[
|\sum_i\log J(w_i;q)-\sum_i\log J(w_i;q_0)|
\le1280Ldh=1280\varepsilon_0=5/4096<1/8.
\]

Consequently, for every w∈A and every matched source radius,

\[
1/2\le\prod_iJ(w_i;q)\le2. \tag{11}
\]

This is a uniform matching statement, not an inference from finitely many q samples.

The full inverse-weighted source shell coefficient is exactly

\[
U(w,\rho)=\frac1{2\ell_*}
\left(\prod_jJ(w_j;q)\right)
\sum_i\frac{N_i}{J(w_i;q)},
\quad\ell_*=1-g/(2R),
\quad N_i=1+\mathbf1_{\{|w_i|<\mathfrak a\}}.
\]

Using 1≤Ni≤2, 1/2≤J≤2, 3/4≤ell*≤1 and (11),

\[
n/8\le U(w,\rho)\le8n. \tag{12}
\]

## 6. Soft-cap loss is paid in the inverse row probability

Normalize each captured source face by its inverse-row density, proportional to 1_capture/G. Its free-coordinate laws have densities proportional to 1_capture/ell_R, with normalizers J≥1/2. They are independent; the pinned source coordinate contributes the fixed log-loss −logell*≤d.

Let V(z)=−logG(z). Since ell≥3/4 and −logell≤(4/3)(1−ell),

\[
\mathbb E_{\rm inv}[-\log\ell_R(Z_i)]
\le2\mathbb E_{\rm full}\frac{-\log\ell_R}{\ell_R}
\le\frac{8g^2}{9R\rho}\le\frac{d^2}{\rho}.
\]

Thus on each captured source face, and therefore their inverse-row mixture,

\[
\mathbb E_{\rm inv}V\le d+nd^2/\rho
\le d+2Ld^2\le3Ld^2.
\]

Markov in **this inverse-row probability**, not the raw source law, gives

\[
\mathbb P_{\rm inv}\{G<\theta\}
=\mathbb P_{\rm inv}\{V>12Ld^2\}\le1/4. \tag{13}
\]

On G≥θ the soft weight b equals θ/G, while everywhere b≤θ/G. Combining (12)–(13), the full soft-weighted shell coefficient lies between

\[
3n\theta/32\quad\text{and}\quad8n\theta. \tag{14}
\]

This loses most of the old raw energy, but recomputes its entire row scale to match the new column scale. There is no exp(Ld³) divisor applied to an unchanged old row budget.

## 7. Cap removal, actual history energy, and the finite-a bridge

Define the actual new kernel as capture, matching, A, source band∩D, multiplied by b(z). Every matching source interval has log-width at least h and at most 4h. The full rows from (14) therefore lie between 3k/32 and 32k.

Because b≤1, removing the original source corner-captured part costs no more than its unweighted full measure. Its mass is at most O(nd exp[−Ld/32]). For d≤2^-24,

\[
\frac{nd e^{-Ld/32}}{nh\theta}
=\frac dh\exp[-Ld/32+12Ld^2]\longrightarrow0.
\]

Hence the final ambient degrees satisfy, eventually,

\[
R_\lambda(w)\ge3k/64\quad(w\in A),
\quad R_\lambda\le32k,
\quad C_\lambda\le4k. \tag{15}
\]

The column bound is (6) and needs no source concentration. Choose any τ≤3k/128. With (10), the same-history energy deficit lemma gives

\[
\Delta\ge\frac{(3k/128)a_*}{32k}
\ge\frac{3nd}{65536}. \tag{16}
\]

For all original finite instances use a≥8n^4, α=a/(a+1), c=(aβ−1)/(a+1), actual μa and its corner envelope λa. On target and source sets, λa=α^-nλ exactly. All actual target items meeting A have radius<3; all source items have radius>4. The guard κ+2/(a+1)<αβ makes every included source precede and capture the same actual target center. Since 0≤b≤1, its weighted actual row is at most τ, with no factor θ^-1. These statements hold μa-a.e. without any distribution hypothesis on μa.

The α^-n scaling multiplies ambient energy by α^-2n and degrees and target mass by α^-n. Since α^-n≥1, (16) remains a lower bound for the finite-a total missing mass, uniformly for all a≥8n^4. It follows that the actual uncaptured load satisfies

\[
U\le\alpha^{-n}J_a(t)-nd/65536. \tag{17}
\]

The weaker rational coefficient in (17) leaves fixed margin below (16).

## 8. Complete scalar cost and the general-L1/all-scales operator conclusion

Put δ=nd/65536 and τ=√[A0(B0−δ)]. The original cap estimates with x∈[L/2,L] give B0>8n eventually, so B0>δ. Also

\[
A_0\le e^x=ne^{-32d^2L},\quad B_0\le nL,
\quad\tau\le n\sqrt L e^{-16d^2L}.
\]

Therefore

\[
\frac\tau k\le\frac{\sqrt L}{h}e^{-4d^2L}\longrightarrow0
\]

at every fixed d>0. In particular τ≤3k/128 eventually. The exact same original covering/selection/witness argument, now with (17), gives for all a≥8n^4

\[
C_n\le\left(1+\frac{(\beta+1/a)^n}{\tau}\right)
\left[\tau+a^{-n}+\alpha^{-n}J_a(t)-\delta\right].
\]

Take a→∞ only in these uniform scalar bounds to obtain (1). Every input is still a general L1 input; no bounded-scale operator, special witness, or new selection history is substituted. The new source weights are used solely in the post-selection energy argument and remain bounded by 1, so all original eligibility and all-positive-common-radius hypotheses are preserved.

For the complete cost, with T=nL−nlogL,

\[
B_0\le1+n\log t
\le T+64nd^2+1,
\]

using −log(1−32d²)≤64d². Further A0+2√(A0B0)=o(n) at each fixed d. For d=2^-30,

\[
\delta=2^{-46}n,\qquad64nd^2=2^{-54}n.
\]

All remaining positive terms are eventually smaller than the fixed margin δ/2−2^-54n. Hence F(A0,B0−δ)≤T−δ/2, proving the proposed (2).

The old and new cost comparison is complete at the same new t: F(A0,B0)−F(A0,B0−δ)>δ. Against older theorems at different t, only eventual asymptotic superiority is claimed; absent explicit thresholds, take the minimum of all applicable bounds.

## 9. Falsifiable sufficient conditions and diagnostic failure modes

The positive mechanism uses four explicit, independently checkable facts: typical matching inverse products (11), the inverse-row log expectation (13), cap loss o(k), and τ/k→0 with the full scalar buffer cost. If any fails, a model reduction in the column maximum alone is not an operator theorem.

Two important invalid shortcuts are excluded. First, applying (13) under the raw source probability would not control row loss: the correct row law is capture/G. Second, choosing x=L−3logL while shrinking rows to θ would generally leave τ/k too large at fixed d. The positive root cost must be recomputed at x=(1−32d²)L; its O(nd²) budget loss is paid explicitly by the larger O(nd) energy deficit.

No explicit n threshold, improved leading nlogn coefficient, or globally optimal d is claimed. The proposed gain is a negative absolute Θ(n) correction, contingent on independent acceptance of the full chain above.

