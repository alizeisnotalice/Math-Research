# R040 draft · Integrated antidiagonal energy and asymptotically sharp degree factors

Status: proposed proof for independent review. It refines the accepted R035 genuine operator mechanism, rather than a finite-instance model cost. The present result is for a common sufficiently-large-n threshold depending on fixed proof parameters. No explicit loglogn≥256 threshold is asserted. No n=7 work, production computation, MATLAB, or Git operation.

## 1. Statement

Write L=logn, l=logL, x=L−3l, t=n/x, β=1+1/t, and let A0=β^n, B0=B_n(t) be the accepted scalar budgets. Fix

\[
0<c<\sqrt2,\qquad0<\epsilon<1/2,
\quad d=c\sqrt{l/L},\quad h=\frac1{Ldl},
\quad I=[\epsilon d,(1-\epsilon)d],
\quad D=(1-2\epsilon)d,
\quad p=L^{-c^2/4},\quad k=nhp.
\]

Let ζ=n^-3 and κ=β−ζ. The letter D in this file is a radial interval length, not the corner-cap domain.

**Proposed theorem.** For any fixed γ satisfying 0<γ<c(1−2ε), and all sufficiently large n,

\[
C_n\le F\left(A_0,B_0-\gamma n\sqrt{l/L}\right)
<W_n-\gamma n\sqrt{l/L}, \tag{1}
\]

where F(A,B)=(√A+√B)^2 and the accepted Wn is

\[
W_n=T+(3l+1)n/L+(9l^2+1)n/L^2,
\quad T=nL-nl.
\]

Consequently for every fixed 0<γ<√2 there is a common sufficiently-large-n threshold such that

\[
C_n\le T-\gamma n\sqrt{l/L}. \tag{2}
\]

To obtain (2) at the stated γ, first use (1) at a slightly larger constant below √2, then absorb Wn−T=o(n√(l/L)). The endpoint γ=√2, a globally optimal proof constant, and any explicit numerical threshold are not claimed.

For the simple fixed choice c=1, ε=1/8, the energy deficit is (3/4−o(1))n√(l/L), so in particular δ=(1/2)n√(l/L) is valid eventually. A finite threshold for that specialization must be certified separately.

## 2. Same ambient geometry, exact kernel and good sets

Use the accepted ambient measure

\[
d\lambda(z)=2^{-n}m(z)^{-n}\mathbf1_{\mathcal D}(z)\,dz,
\quad m(z)=\max(1,\|z\|_\infty),
\]

\[
\mathcal D=\{\|z\|_\infty<t:\exists i,\ z_i<1-\beta m(z)\}.
\]

Before cap restriction, each face at radius R>1 has full-shell measure dR/(2R), with n−1 independent Uniform[−R,R] free coordinates. Full shell mass is n dR/R.

Targets have r=||w||∞=1+u, u∈I. Sources have ρ=||z||∞=tq, q=1−v, v∈I. Matching is |u+v−d|<h. For all sufficiently large n, r<2, ρ>4, q>1/2, 1<κ<β, and h+ζt<d/2.

For targets set

\[
S_w=\sum_i(|w_i|-(r-d))_+,
\quad m_w=d+(n-1)d^2/(2r),
\quad K_w=\#\{i:|w_i|>r-2d\},
\]

\[
\mathcal A=\mathcal D\cap\{u\in I,
|S_w-m_w|\le t/l,\ K_w\le4nd\}.
\]

For sources set

\[
S_z=\sum_i(|z_i|-(\rho-d))_+,
\quad m_z=d+(n-1)d^2/(2\rho),
\quad K_z=\#\{i:|z_i|>\rho-2d\},
\]

\[
\mathcal P=\mathcal D\cap\{v\in I,
|S_z-m_z|\le1/l,\ K_z\le8Ld\}.
\]

All faces of the target annulus are used; this is not a negative-face-only target. Define the Borel kernel

\[
K(w,z)=\mathbf1_{\mathcal A}(w)\mathbf1_{\mathcal P}(z)
\mathbf1_{\{|u+v-d|<h\}}
\mathbf1_{\{\|z-w\|_\infty<\kappa\rho\}}. \tag{3}
\]

A and P are disjoint radial sets. All degrees below are λ-degrees, never assumed bounds on the unknown actual μ.

## 3. Good probabilities and deletion costs

For a free coordinate X uniform on [−R,R],

\[
\mathbb E(|X|-(R-d))_+=d^2/(2R),\quad
\mathbb E[(|X|-(R-d))_+^2]=d^3/(3R).
\]

On a fixed face the exceptional coordinate contributes d deterministically. Thus the displayed mw,mz are the exact face-conditional means.

For targets, uniformly in r and in the chosen maximal face,

\[
\operatorname{Var}S_w\le nd^3/(3r),
\quad\mathbb P(|S_w-m_w|>t/l)
=O(L^2d^3l^2/n)=o(1).
\]

Furthermore Kw=1+Bin(n−1,2d/r), its mean is (2+o(1))nd and its variance O(nd). Chebyshev at 4nd gives bad probability O(1/(nd))=o(1), improving the fixed Markov loss in R035. Hence target good probability is 1−o(1), uniformly on every face and radius. The corner-captured target fraction is exponentially small in n, because ell_t(r)≤3/4 eventually. Consequently the target radial density satisfies

\[
d\lambda_{\mathcal A}(u)=n[1-o(1)]\frac{du}{1+u}
\quad\text{uniformly on }I,
\qquad a_*:=\lambda(\mathcal A)=nD[1+o(1)]. \tag{4}
\]

For sources,

\[
\operatorname{Var}S_z\le nd^3/(3\rho)=O(Ld^3),
\quad\mathbb P(|S_z-m_z|>1/l)=O(Ld^3l^2).
\]

Kz=1+Bin(n−1,2d/ρ), with mean (2+o(1))Ld and variance O(Ld), so the bad probability at 8Ld is O(1/(Ld)). Multiplying by full source annulus mass O(nd), the total bad-statistics mass is

\[
b_{\rm stat}=O(nLd^4l^2+n/L)
=O_c(nl^4/L+n/L). \tag{5}
\]

The source corner-captured part has coefficient at most (n/2)exp(−aεLd) for an absolute a>0, because v≥εd and x/L→1. Its mass is O(nd exp(−aεLd)). Therefore the total source mass removed from a full annulus satisfies

\[
b_n=O_c(nl^4/L+n/L)+O(nd e^{-a\epsilon Ld})=o(k). \tag{6}
\]

For clarity, the leading ratio in (6) is

\[
O_c\left(l^{11/2}L^{-1/2+c^2/4}
+l^{3/2}L^{-1/2+c^2/4}\right)\longrightarrow0,
\]

using c²<2. Deleting b_n in a row costs at most b_n because K≤1. Later, integrated energy loss is at most λ(A)b_n=o(nDk), so no unjustified conditional independence or bad-column degree estimate is required.

## 4. Relative Gaussian product asymptotics

For a matched pair of radial coordinates put e=u+v−d=r−q−d, |e|<h, and b=e+ζρ. The threshold displacement has size at most h+ζt=o(d).

For w∈A, normalized full-source coordinate lengths are

\[
L_i^{\rm row}=1-
\frac{(|w_i|-(r-d-b))_+}{2\rho}.
\]

Only coordinates counted by Kw can change when moving the threshold from r−d by b. Consequently

\[
\left|\sum_i(1-L_i^{\rm row})-
\frac{nd^2}{4r\rho}\right|
\le\frac{t/l+4nd(h+\zeta t)+2d}{2\rho}
=O(1/l)+o(1)=o(1). \tag{7}
\]

For z∈P, normalized full-target lengths are

\[
L_i^{\rm col}=1-
\frac{(|z_i|-(\rho-d-b))_+}{2r},
\]

and Kz gives

\[
\left|\sum_i(1-L_i^{\rm col})-
\frac{nd^2}{4\rho r}\right|
\le\frac{1/l+8Ld(h+\zeta t)+2d}{2r}
=O(1/l)+o(1)=o(1). \tag{8}
\]

Uniformly across the radial strip,

\[
\frac{nd^2}{4r\rho}
=\frac{c^2l}{4}\frac{x}{Lrq}
=c^2l/4+O_c(ld+l^2/L)=c^2l/4+o(1). \tag{9}
\]

The maximum row loss is O(d/t), the maximum column loss O(d), and each total loss is O(l). Thus the sums of squared losses are o(1). Removing one coordinate also changes the logarithm by o(1). The elementary log(1−s) bounds now imply, uniformly over all good points, all matched radii and every deleted coordinate j,

\[
\prod_{i\ne j}L_i^{\rm row}
=p[1+o(1)],\qquad
\prod_{i\ne j}L_i^{\rm col}
=p[1+o(1)]. \tag{10}
\]

These are relative asymptotics, improving the fixed e^±2 enclosure in R035. The o(1) terms are absolute once c,ε are fixed.

## 5. Face factors approach n, not n/2

The full-source face coefficient at a fixed good w is exactly

\[
\frac12\sum_i
\left[1+\mathbf1_{\{|w_i|<(\kappa-1)\rho\}}\right]
\prod_{j\ne i}L_j^{\rm row}.
\]

Any missing second face has |wi|≥r−d−b>r−2d, so their count is at most Kw≤4nd=o(n). With (10), the whole coefficient is

\[
np[1+o(1)] \tag{11}
\]

uniformly for all matched radii and w∈A.

For fixed z∈P the full-target coefficient is

\[
\frac12\sum_i
\left[1+\mathbf1_{\{|z_i|<\kappa\rho-r\}}\right]
\prod_{j\ne i}L_j^{\rm col}.
\]

Missing second faces are contained in Kz≤8Ld=o(n), so this coefficient also equals np[1+o(1)]. No target-negative-face factor 1/2 is lost, since the full target annulus is used. Restricting the target to A only decreases a column upper bound.

## 6. Integrated strip energy removes the endpoint factor

Let Rλ(w), Cλ(z) be the row and column degrees of (3). Each radial matching interval has length at most 2h, so (11) and its column counterpart give

\[
M_n:=\max\{\sup_{w\in A}R_\lambda(w),
\sup_{z\in P}C_\lambda(z)\}
\le2nhp[1+o(1)]=2k[1+o(1)]. \tag{12}
\]

For the energy lower bound, do not use the shortest endpoint interval h. Instead integrate over the whole square I×I. Its antidiagonal strip has exact Euclidean area

\[
|\{(u,v)\in I^2:|u+v-d|<h\}|=2hD-h^2, \tag{13}
\]

for all sufficiently large n when h<D. The Jacobian 1/[(1+u)(1−v)] equals 1+O(d), uniformly. From target radial density (4) and full-source coefficient (11), energy before source deletion is

\[
E_0=n^2p(2hD-h^2)[1+o(1)].
\]

Actual source-good and cap deletion costs at most a_*b_n=o(nDk) by (6). Since h/D=1/[(1−2ε)c²l²]→0, the final ambient energy obeys

\[
E_{\lambda\lambda}=\iint Kd\lambda d\lambda
\ge2n^2hpD[1-o(1)]=2knD[1-o(1)]. \tag{14}
\]

Thus endpoints do not cost a permanent factor 1/2. Both the energy and maximum degrees use the same kernel and the same full ambient measure.

## 7. Energy deficit and exact original-history transfer

Use the already accepted bipartite deficit lemma: for μ≤λ with A/P disjoint, actual rowμ≤τ, and ambient row/column degrees≤M,

\[
\Delta=\lambda(H)-\mu(H)
\ge\frac{E_{\lambda\lambda}-\tau\lambda(A)}{M}.
\]

For any τ≤n/L,

\[
\frac\tau k\le c l^{3/2}L^{-1/2+c^2/4}\longrightarrow0,
\]

so (4), (12) and (14) imply

\[
\Delta\ge nD[1-o(1)]. \tag{15}
\]

For completeness, this remains an estimate for all actual finite histories, not just a continuum model. In the original normalization, μa is the actual uncaptured-increment measure and λa the accepted corner envelope. On A and P, λa=α^-nλ exactly because D∞⊂Ha. Choose a≥8n^4. Then κ+2/(a+1)<αβ. Every actual target whose increment lies in A has radius <3; every source from P has radius >4, so all sources precede that same target and their actual centers are captured. Thus rowμa≤τ holds μa-a.e. on A, even after radial matching restrictions.

If ca=α^-n≥1, ambient degrees scale by ca, energy by ca² and target mass by ca. The finite-a deficit is at least

\[
\frac{c_a^2E_{\lambda\lambda}-\tau c_aa_*}{c_aM_n}
\ge\frac{E_{\lambda\lambda}-\tau a_*}{M_n}
\ge nD[1-o(1)],
\]

uniformly for every real a≥8n^4. The unknown μa is not assumed non-concentrated. No selector limit occurs.

Fix γ<c(1−2ε). Since nD=c(1−2ε)n√(l/L), (15) gives deficit≥δ=γn√(l/L) for all sufficiently large n. Let τ=√[A0(B0−δ)]. The accepted estimates B0>8n, B0≤nL and A0≤n/L³ give δ<B0 and τ≤n/L eventually. The original finite cover/load formula, with this deficit, then yields

\[
C_n\le\left(1+\frac{(\beta+1/a)^n}{\tau}\right)
\left[\tau+a^{-n}+\alpha^{-n}J_a(t)-\delta\right].
\]

Taking a→∞ only in these uniform scalar inequalities proves (1). It preserves the original general-L1 input and all-positive-common-radius scope. The accepted Wn bound and the root difference >δ give the second inequality in (1).

Finally, for any γ<√2 choose c∈(γ,√2) and ε small enough that c(1−2ε)>γ. This proves the operator mechanism at every fixed coefficient below √2. Use an intermediate coefficient above γ and Wn−T=o(n√(l/L)) to conclude (2).

## 8. What remains quantitative

The theorem above has a common sufficiently-large-n threshold depending on fixed c,ε,γ. The elementary rates in (5)–(10) provide a direct route to quantification, but no concrete threshold has been certified in this draft. In particular the suggested c=1,ε=1/8,γ=1/2 statement at loglogn≥256 remains a separate task.

This is not a proof that √2 is the optimal constant for all covering arguments. It is the endpoint approached by this specific Gaussian degree calculation while τ/k→0 and bad-source/k→0. It does not lower the leading coefficient 1 of nlogn. No finite matrix, optimizer, or unquantified numerical threshold replaces the uniform analytic argument.

