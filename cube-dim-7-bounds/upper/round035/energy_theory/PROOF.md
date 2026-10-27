# R035 draft · Antidiagonal multiscale energy cut

Status: complete proposed proof, submitted for independent audit; not yet a frozen result. All asymptotic statements below have absolute constants uniform over the displayed radial intervals and good sets. The theorem claims a common sufficiently-large-n threshold, not a presently certified numerical threshold. No n=7 work, old grid replay, MATLAB, or Git operations.

## 1. Proposed theorem and exact parameters

Let n be an integer, L=log n, l=log L, x=L−3l, t=n/x, β=1+1/t. Write A0=β^n and B0=B_n(t) for the accepted original scalar budget. Define

\[
d=\sqrt{l/L},\qquad h=\frac1{100Ld},\qquad
I=[d/4,3d/4],\qquad
\zeta=n^{-3},\quad\kappa=\beta-\zeta,
\quad k=n h L^{-1/4}.
\]

For all sufficiently large n, the proposed conclusion is

\[
C_n\le F\left(A_0,B_0-\frac{nd}{65536}\right)
<W_n-\frac{nd}{65536}, \tag{1}
\]

where F(A,B)=(√A+√B)^2 and the accepted W_n equals

\[
T+(3l+1)n/L+(9l^2+1)n/L^2,\qquad T=nL-nl.
\]

Consequently, for all sufficiently large n,

\[
C_n\le n\log n-n\log\log n
-\frac{n}{131072}\sqrt{\frac{\log\log n}{\log n}}. \tag{2}
\]

This would change the order of the R033/034 negative correction, while retaining leading coefficient 1.

## 2. Ambient measures and exact radial coordinates

Put m(z)=max(1,||z||∞), and let

\[
D_\infty=\{\|z\|_\infty<t:\exists i,\ z_i<1-\beta m(z)\},
\qquad
d\lambda(z)=2^{-n}m(z)^{-n}\mathbf1_{D_\infty}(z)\,dz.
\]

For estimates before removing the corner-captured part, use the full-shell measure

\[
d\lambda^0(z)=2^{-n}\|z\|_\infty^{-n}\mathbf1_{\{\|z\|_\infty>1\}}\,dz.
\]

At radius R>1, every face z_i=±R has measure dR/(2R), and its remaining n−1 coordinates are independent Uniform[−R,R]. Thus the full shell has measure n dR/R. Ties of maximal coordinates are null.

The target annulus has radius r=1+u, u∈I. The source annulus has radius ρ=tq, q=1−v, v∈I. They are disjoint, with r<2 and ρ>4 for large n. Matching is antidiagonal:

\[
|u+v-d|<h,\quad\text{equivalently }|r-q-d|<h. \tag{3}
\]

For fixed u∈I the allowed v interval has length between h and 2h; likewise after fixing v. This follows because its center d−u belongs to I. Also dρ/ρ=dv/(1−v) in absolute value, between dv and 2dv, while dr/r=du/(1+u)≤du. The full target shell is used, not only its negative faces.

## 3. Good sets

For a target point w with r=||w||∞=1+u, set

\[
S_w=\sum_i(|w_i|-(r-d))_+,
\quad m_w=d+(n-1)d^2/(2r),
\quad K_w=\#\{i:|w_i|>r-2d\}.
\]

Define

\[
\mathcal A= D_\infty\cap\{u\in I,
|S_w-m_w|\le t/100,\ K_w\le12nd\}.
\]

For a source z with ρ=||z||∞=t(1−v), set

\[
S_z=\sum_i(|z_i|-(\rho-d))_+,
\quad m_z=d+(n-1)d^2/(2\rho),
\quad K_z=\#\{i:|z_i|>\rho-2d\}.
\]

Define

\[
\mathcal P=D_\infty\cap\{v\in I,
|S_z-m_z|\le1/100,\ K_z\le16Ld\}.
\]

All sets are Borel. Define the exact kernel

\[
K(w,z)=\mathbf1_{\mathcal A}(w)\mathbf1_{\mathcal P}(z)
\mathbf1_{\{|r-q-d|<h\}}
\mathbf1_{\{\|z-w\|_\infty<\kappa\rho\}}. \tag{4}
\]

All row and column degrees in this proof are ambient λ-degrees of this kernel; they are not μ-degrees.

## 4. Good-set mass estimates

For a free coordinate X uniform on [−R,R] and 0<d<R,

\[
\mathbb E(|X|-(R-d))_+=d^2/(2R),\qquad
\mathbb E[(|X|-(R-d))_+^2]=d^3/(3R).
\]

On any fixed face one coordinate contributes d deterministically. Hence the centers m_w,m_z above are exact conditional means.

For target faces, Var S_w≤nd^3/(3r), and Chebyshev bounds its bad probability by O(L²d³/n)=o(1). Also E K_w≤1+2nd≤3nd, so Markov gives P(K_w>12nd)≤1/4. Therefore each target face has at least half good points for large n. The target annulus has full mass

\[
n\log\frac{1+3d/4}{1+d/4}\ge nd/4.
\]

The corner-captured target part is at most O(nd(3/4)^{n-1}): its shell coefficient is (n/2)ell_t(r)^{n−1}, and ell_t(r)≤3/4 for large n. Thus

\[
a_*=\lambda(\mathcal A)\ge nd/16. \tag{5}
\]

For source faces, q≥1/2 and ρ>2d give

\[
\operatorname{Var}S_z\le nd^3/(3\rho)\le(2/3)Ld^3.
\]

The bad S_z probability is at most (20000/3)Ld³. Moreover K_z=1+Bin(n−1,2d/ρ) on each face. For Ld≥1,

\[
\mathbb E K_z\le5Ld,\quad\operatorname{Var}K_z\le4Ld,
\quad\mathbb P(K_z>16Ld)\le4/(121Ld).
\]

The full source annulus has mass at most nd, so its bad-statistics mass is

\[
\le\frac{20000}{3}\frac{nl^2}{L}+\frac4{121}\frac nL=o(k). \tag{6}
\]

Indeed k=n/[100 L^(3/4)√l], so the ratio tends to zero like O(l^(5/2)L^(-1/4)). The corner-captured source part is also o(k): throughout v∈I, its coefficient is at most

\[
\frac n2\exp(-Ld/32),
\]

because (n−1)(1/ρ−1/t)/2≥Ld/32. Multiplying by the band log-length≤d gives O(nd exp(−Ld/32))=o(k). Consequently removing both the source bad sets and its corner-captured part costs o(k) in any row, even if removal is maximally concentrated. No conditional independence of good sets is assumed.

## 5. Uniform Gaussian product bounds

This section concerns any matched r,q satisfying (3). Put e=r−q−d, |e|<h, and b=e+ζρ. Then |b|≤h+ζt=o(d). The following estimates are uniform over all good points, all matched radii, and all deleted-coordinate products.

### 5.1 Row products

For a fixed w∈A, the full source capture box has normalized coordinate lengths

\[
L_i^{\rm row}=\min\left(1,\frac{1+\kappa}{2}
-\frac{|w_i|}{2\rho}\right)
=1-\frac{(|w_i|-(r-d-b))_+}{2\rho}.
\]

Since |b|<d, changing the threshold from r−d only affects coordinates counted by K_w. Thus the capture exceedance sum differs from S_w by at most 12nd(h+ζt). If T_w=Σi(1−L_irow),

\[
\left|T_w-\frac{nd^2}{4r\rho}\right|
\le\frac{t/100+12nd(h+\zeta t)+2d}{2\rho}.
\]

The limsup of the right side is at most 1/200+3/50=13/200. Uniformly,

\[
\frac{nd^2}{4r\rho}
=\frac l4\frac{x}{Lrq}=l/4+o(1),
\]

because l d+l²/L→0. Also max_i(1−L_irow)=O(d/t), so Σi(1−L_irow)^2=o(1). The elementary bound −v−v²/(1−v)≤log(1−v)≤−v then gives, for all sufficiently large n and every deleted coordinate j,

\[
e^{-2}L^{-1/4}\le\prod_{i\ne j}L_i^{\rm row}
\le e^2L^{-1/4}. \tag{7}
\]

### 5.2 Column products

For a fixed z∈P and matched target radius r, the full target capture box has lengths

\[
L_i^{\rm col}=\min\left(1,\frac{r+\kappa\rho-|z_i|}{2r}\right)
=1-\frac{(|z_i|-(\rho-d-b))_+}{2r}.
\]

Threshold movement changes the exceedance sum by at most 16Ld(h+ζt). Therefore, for T_z=Σi(1−L_icol),

\[
\left|T_z-\frac{nd^2}{4\rho r}\right|
\le\frac{1/100+16Ld(h+\zeta t)+2d}{2r}.
\]

Its uniform limsup is at most (1/100+4/25)/2=17/200. The same Gaussian mean is l/4+o(1); now max loss is O(d), and Σ loss²=O(dl)=o(1). Thus for every deleted coordinate j,

\[
e^{-2}L^{-1/4}\le\prod_{i\ne j}L_i^{\rm col}
\le e^2L^{-1/4} \tag{8}
\]

for all sufficiently large n. Equations (7)–(8) use the exact centering under fixed-face coordinates; no unconditioned cube law replaces a face law.

## 6. Ambient row and column degrees

Let R(w)=∫K(w,z)dλ(z), C(z)=∫K(w,z)dλ(w). For κ>1 and ρ>s, a source coordinate always has at least one face captured, and at most two. The full source shell coefficient is

\[
\frac12\sum_i
\left[1+\mathbf1_{\{|w_i|<(\kappa-1)\rho\}}\right]
\prod_{j\ne i}L_j^{\rm row}.
\]

By (7), it is between (n/2)e^-2 L^-1/4 and n e² L^-1/4. Matching gives log-length between h and 4h. Removing the o(k) source mass from Section 4 leaves, uniformly on A,

\[
r_*:=e^{-2}k/4\le R(w)\le4e^2k. \tag{9}
\]

For a fixed z, the full target shell coefficient is

\[
\frac12\sum_i
\left[1+\mathbf1_{\{|z_i|<\kappa\rho-r\}}\right]
\prod_{j\ne i}L_j^{\rm col}.
\]

It is at most n e² L^-1/4 by (8). Matching target log-length is at most 2h, and restricting to target good sets and cap only decreases mass. Hence

\[
C(z)\le2e^2 k. \tag{10}
\]

Notice that full target faces are used. If only negative target faces were used, this formula would need their one-sided indicators and an additional factor 1/2 in suitable comparisons.

## 7. Energy deficit lemma and the counterexample to free addition

Let λ be any finite ambient measure, μ≤λ, and let A,P be disjoint. Suppose 0≤K≤1 is supported in A×P and actual history yields

\[
\int K(w,z)d\mu(z)\le\tau\quad\text{for }\mu\text{-a.e. }w\in A.
\]

Write Eνν=∫∫K dνdν, d0=λ−μ, Δ=d0(total). If Rλ≤M and Cλ≤M, then

\[
\Delta\ge\frac{E_{\lambda\lambda}-\tau\lambda(A)}{M}. \tag{11}
\]

Indeed Eμμ≤τμ(A)≤τλ(A), while

\[
E_{\lambda\lambda}-E_{\mu\mu}
=E_{d_0|A,\lambda|P}+E_{\mu|A,d_0|P}
\le M[d_0(A)+d_0(P)]\le M\Delta.
\]

Thus μ need not be spread out; degrees belong to λ. This is a genuine common-history estimate, not a sum of independent worst positions.

Counterexample diagnostic: take disjoint Lebesgue intervals A of mass N and P of mass 1, K=1 on A×P, μ=λ on A and μ=τλ on P, 0<τ<1. Every target's local inequality is saturated, but the total deficit is just 1−τ, however large N is. Independently adding N cuts would be false. Formula (11) is sharp here because the column degree equals N. This is an abstract shared-capacity counterexample, not a claimed cube history or a lower bound for Cn.

Apply (11) to the degrees (9)–(10), with M=4e²k. Since

\[
\frac{k}{n/L}=\frac{L^{1/4}}{100\sqrt l}\longrightarrow\infty,
\]

we have τ≤n/L≤r*/2 for large n. Therefore (5) and (9) give

\[
\Delta\ge\frac{a_*r_*}{2M}
\ge\frac{e^{-4}}{512}nd>\frac{nd}{65536}. \tag{12}
\]

The final rational comparison follows from e<3, hence e^4<81<128.

## 8. Exact finite-a original-history bridge

Use the accepted finite-a actual increment measure μa and ambient λa from R033:

\[
d\lambda_a=\alpha^{-n}2^{-n}m^{-n}\mathbf1_{H_a}dz,
\quad H_a=\{\|z\|<t:\exists i,z_i\le\alpha-cm(z)\},
\quad\alpha=a/(a+1),\ c=(a\beta-1)/(a+1).
\]

For every z∈D∞, α−cm−(1−βm)=((β+1)m−1)/(a+1)>0. Consequently λa=α^-nλ exactly on A and P. Actual uncaptured increment measure satisfies μa≤λa; no new density hypothesis is imposed on μa.

Take a≥8n^4. Then a>2t, κ>1, and κ+2/(a+1)<αβ. If w∈A belongs to actual item j, its radius satisfies rj≤||u−w||/c<α+s<3, where s=1+3d/4<2. Any matched source z∈P belongs to an item i of radius ri≥ρ>4; therefore i<j. Center displacements imply capture by bQi, since

\[
\|v_i-v_j\|<\|z-w\|+(r_i+r_j)/(a+1)
< [\kappa+2/(a+1)]r_i<\alpha\beta r_i.
\]

Thus for μa-a.e. w∈A, its one actual eligibility pays ∫K(w,z)dμa(z)≤τ. Matching merely restricts sources, and no sources from other histories are introduced.

On A/P the ambient scaling factor is c_a=α^-n≥1. Their row and column upper degrees are bounded by c_aM, their energy is c_a²Eλλ, and their target mass is c_aa*. Applying (11) gives a total original-envelope deficit at least

\[
\frac{c_a^2E_{\lambda\lambda}-\tau c_aa_*}{c_aM}
\ge\frac{E_{\lambda\lambda}-\tau a_*}{M}
>nd/65536.
\]

This holds for every real a≥8n^4 and every original finite history, not just along a convergent family.

Set δ=nd/65536 and τ=√[A0(B0−δ)]. The accepted old bounds B0>8n, B0≤nL, A0≤n/L³ guarantee δ<B0 and τ≤n/L for sufficiently large n. The finite cover/load formula with the just-proved deficit yields

\[
C_n\le\left(1+\frac{(\beta+1/a)^n}{\tau}\right)
\left[\tau+a^{-n}+\alpha^{-n}J_a(t)-\delta\right].
\]

Let a→∞ only in these uniform scalar inequalities. This gives (1); subtraction of the accepted Wn bound gives its strict comparison. Finally its positive remainder is O(nl/L)=o(nd), proving (2). Neither a selector limit nor a global parameter minimizer is needed.

## 9. Remaining audit status

The proposed theorem changes order only after all uniform product, bad-mass, and finite-a estimates above are accepted. No numerical matrix or asymptotic heuristic is substituted for those estimates. Independent auditors have been given this exact path. A numerical threshold log L≥512 is plausible from the displayed rates but is not asserted here until separately certified; the theorem as drafted requires only a common sufficiently-large-n threshold.

