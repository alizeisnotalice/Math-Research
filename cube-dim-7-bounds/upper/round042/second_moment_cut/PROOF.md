# R042 draft · Second-moment column truncation at a new scalar scale

Status: proposed full operator proof, submitted for independent review. It is independent of the fixed-x=L−3logL R040 refinement, whose source remains frozen. The threshold here is a common sufficiently-large-n threshold, not a certified finite numerical value. No n=7 work, production computation, MATLAB, or Git operation.

## 1. Proposed conclusion and fixed rational constants

Let

\[
L=\log n,\quad s=L^{1/3},\quad x=L-s,\quad t=n/x,
\quad\beta=1+1/t,
\quad d=L^{-1/3}=1/s,\quad h=\frac1{100Ld}=\frac{d^2}{100},
\]

\[
I=[d/4,3d/4],\quad p=e^{-s/4},\quad M_0=nhp,
\quad\zeta=n^{-3},\quad\kappa=\beta-\zeta.
\]

Use the original scalar quantities A0=β^n and B0=B_n(t). Fix

\[
c_E=2^{-10},\quad C_R=36,\quad C_2=144,\quad H=2^{20}.
\]

The proposed result, for every sufficiently large n, is

\[
C_n\le F\left(A_0,B_0-2^{-31}nd\right), \tag{1}
\]

and consequently

\[
C_n\le n\log n-n\log\log n
-2^{-32}\frac{n}{(\log n)^{1/3}}. \tag{2}
\]

The constants are intentionally conservative exact rational numbers. No claim is made that this bound dominates every previously applicable finite-dimensional bound; the valid reported bound is their minimum. Its asymptotic correction has a larger order than n√(loglogn/logn).

## 2. Ambient measure, target set and raw kernel

Let

\[
\mathcal D=\{\|z\|_\infty<t:\exists i,z_i<1-\beta m(z)\},
\quad m(z)=\max(1,\|z\|_\infty),
\quad d\lambda=2^{-n}m(z)^{-n}\mathbf1_{\mathcal D}(z)dz.
\]

The full-shell measure before cap restriction has each face at radius R weighted by dR/(2R), with its other n−1 coordinates independent Uniform[−R,R]. Thus full shell mass is n dR/R.

Targets have radius r=||w||∞=1+u, u∈I. Sources have radius ρ=||z||∞=tq, q=1−v, v∈I. Matching means

\[
|u+v-d|<h,\qquad e=r-q-d\in(-h,h). \tag{3}
\]

For targets put

\[
S_w=\sum_i(|w_i|-(r-d))_+,
\quad m_w=d+(n-1)d^2/(2r),
\quad K_w=\#\{i:|w_i|>r-2d\},
\]

\[
\mathcal A=\mathcal D\cap\{u\in I,
|S_w-m_w|\le t/100,\ K_w\le12nd\}.
\]

Let Praw= D∩{v∈I}. No source good-statistics condition is imposed. The raw kernel Kraw is the indicator of A×Praw, matching (3), and ||z−w||∞<κρ. We first bound its row and energy; a source degree truncation is introduced later.

All eventual guards used below are elementary: d→0, t∼n/L→∞, r<2, ρ>4, q>1/2, h+ζt<d/2, 1<κ<β, and x∈[L/2,L].

## 3. Target mass and uniformly bounded row-product errors

For X uniform on [−r,r], the excess over r−d has mean d²/(2r) and second moment d³/(3r). On a fixed maximal face the fixed coordinate contributes d. Thus mw is the exact conditional mean and VarSw≤nd³/(3r). Chebyshev at t/100 gives bad probability O(L²d³/n)=O(L/n)=o(1). Also E Kw≤1+2nd≤3nd, so Markov at 12nd loses at most 1/4. Hence every target face retains at least half its points eventually.

The full target annulus mass is at least nd/4. The original corner-captured target part is O(nd(3/4)^(n−1)), because ell_t(r)≤3/4 eventually. Therefore

\[
a_*:=\lambda(\mathcal A)\ge nd/16,
\qquad a_*\le nd. \tag{4}
\]

For a matched source radius define b=e+ζρ. The normalized source coordinate capture lengths at fixed w are

\[
L_i(w,\rho)=1-
\frac{(|w_i|-(r-d-b))_+}{2\rho}.
\]

Since |b|<d, threshold movement is supported on coordinates counted by Kw. The sum T_w=Σi(1−Li) satisfies

\[
\left|T_w-\frac{nd^2}{4r\rho}\right|
\le\frac{t/100+12nd(h+\zeta t)+2d}{2\rho}.
\]

Its uniform limsup is at most 13/200, since x/L→1 and Ld³=1. The Gaussian center is now only within a bounded error, not o(1), of s/4:

\[
\frac{nd^2}{4r\rho}=\frac{s-d}{4rq},
\quad |rq-1|\le d,
\quad\left|\frac{s-d}{4rq}-\frac s4\right|
\le\frac{1+d}{4(1-d)}.
\]

The right side tends to 1/4. Maximum coordinate loss is O(d/t); the sum of squared losses is O(1/t)=o(1), and deleting one coordinate changes the logarithm by o(1). Thus for every omitted coordinate j, uniformly on all good targets and matched radii,

\[
e^{-2}p\le\prod_{i\ne j}L_i(w,\rho)\le e^2p \tag{5}
\]

eventually. The constants e^±2 are deliberately wider than the displayed limiting errors; no bounded error is silently replaced by o(1).

Each source coordinate has at least one captured face and at most two because κ>1 and ρ>r. Thus full-source face coefficients lie between (n/2)e^-2p and ne²p. Every target's matching v interval has length between h and 2h, so its source log-length lies between h and 4h. The original corner-captured source part has total mass

\[
O(nd\exp[-aLd])=o(M_0)
\]

for an absolute a>0; here Ld=L^(2/3), which dominates s. Removing that cap complement costs at most its entire mass in any row. Hence raw degrees satisfy

\[
R_{\rm raw}(w)\ge e^{-2}M_0/4,
\qquad R_{\rm raw}(w)\le4e^2M_0<36M_0. \tag{6}
\]

The raw ambient energy therefore obeys

\[
E_{\rm raw}\ge a_*e^{-2}M_0/4
\ge e^{-2}ndM_0/64>c_EndM_0, \tag{7}
\]

using e<3 and cE=1/1024. No source mass-deletion argument will be used for the subsequent column cutoff.

## 4. A pointwise column envelope without source concentration

For each source z with ρ=tq and v∈I set

\[
r_-=q+d-h>1,
\quad g=r_--q+\zeta\rho=d-h+\zeta\rho.
\]

Every matched target radius r lies at least r_-. This lower radius need not itself lie in the target annulus; that only makes it a conservative bound.

For fixed z and target radius r, put ai=κρ−|zi|>0. The normalized target side length is

\[
\ell_i(r)=\min(1,(1+a_i/r)/2).
\]

It is nonincreasing in r. All these lengths are at least 1/2. A full target-shell capture coefficient is at most 2n times the product of its n side lengths, because each deleted-coordinate product is at most twice the full product, and the normalized face factor is 1/2. The matched target interval has ordinary width≤2h and logarithmic width≤2h/r_-. Hence the full target column, and therefore its restriction to A, obey

\[
C_{\rm raw}(z)\le Y(z):=
\frac{4nh}{r_-}\prod_{i=1}^n\ell_i(r_-). \tag{8}
\]

This argument is pointwise in ρ. It never differentiates r_-(ρ), and interval endpoints can only shrink the actual matched column. A sharper coefficient bound is possible, but is not needed.

## 5. Exact conditional second moment, including the pinned face

For all sufficiently large n, ζt≤h/2, so 0<d−h≤g≤d. At radius ρ the source has one pinned coordinate |zj|=ρ. Its side length is exactly

\[
\ell_*=1-g/(2r_-)\in(1/2,1),
\]

not an independent uniform coordinate. Its squared factor can be discarded only in the upper-bound direction. For every other coordinate X=|Zi| uniform on [0,ρ],

\[
\ell(X)=1-\frac{(X-(\rho-g))_+}{2r_-},
\]

and exact integration gives

\[
\mathbb E\ell(X)^2
=1-\frac{g^2}{2r_-\rho}
+\frac{g^3}{12r_-^2\rho}. \tag{9}
\]

The following coarse absolute bounds suffice. Since r_-q≤1+3d/4, h=d²/100 and xd²=s−d,

\[
\frac{(n-1)g^2}{2r_-\rho}
\ge\frac{s-d}{2}(1-1/n)(1-d/50)(1-3d/4)
\ge s/2-1 \tag{10}
\]

eventually. For example the loss from s/2 is bounded by d/2+s/(2n)+77/200, which is <1 eventually. Also, using r_->1 and q≥1/2,

\[
\frac{(n-1)g^3}{12r_-^2\rho}
\le\frac{xd^3}{12q}\le1/6. \tag{11}
\]

Apply logu≤u−1 to (9). Equations (10)–(11) imply

\[
(\mathbb E\ell^2)^{n-1}
\le e^{-s/2+7/6}<9p^2.
\]

Here e^(7/6)<e²<9. Conditional on every maximal source face,

\[
\mathbb E Y^2
=\left(\frac{4nh}{r_-}\right)^2
\ell_*^2(\mathbb E\ell^2)^{n-1}
\le144M_0^2. \tag{12}
\]

The full source annulus has measure≤nd, so

\[
\int Y^2\,d\lambda_{\rm source}^0
\le C_2ndM_0^2,\qquad C_2=144. \tag{13}
\]

Restriction to the original corner cap only decreases this moment. The O(1) compensation in (10)–(12) is retained; it is not asserted to vanish.

## 6. Truncate columns by energy, not by deleted source mass

Define

\[
\mathcal P=\mathcal P_{\rm raw}\cap\{Y\le HM_0\},
\qquad H=2^{20},
\]

and K as the raw kernel restricted to A×P. Its column degree is at most HM0, and its row degree remains at most 36M0≤HM0.

The removed energy is exactly the integral of raw target columns over the removed sources. By (8) and (13),

\[
E_{\rm raw}-E_{\rm good}
\le\int_{Y>HM_0}Y\,d\lambda
\le\frac{1}{HM_0}\int Y^2d\lambda
\le\frac{C_2}{H}ndM_0
\le\frac{c_E}{4}ndM_0. \tag{14}
\]

The last inequality is the exact integer comparison H≥4C2/cE. This does not require the deleted source mass to be o(M0), and does not assume any concentration of the unknown actual measure μ. Thus

\[
E_{\rm good}\ge(3c_E/4)ndM_0,
\qquad\max\{R_\lambda,C_\lambda\}\le HM_0. \tag{15}
\]

## 7. Energy deficit and all original finite histories

The accepted energy lemma says that if μ≤λ, A/P are disjoint, and actual eligibility gives rowμ≤τ on μ-a.e. target, then

\[
\Delta\ge\frac{E_{\lambda\lambda}-\tau\lambda(A)}
{\max(R_\lambda,C_\lambda)}.
\]

Whenever τ/M0≤cE/4, equations (4) and (15) give

\[
\Delta\ge\frac{(c_E/2)ndM_0}{HM_0}
=2^{-31}nd. \tag{16}
\]

This transfers exactly to the original finite-a algorithm. Let μa be its actual uncaptured-increment measure and λa its corner envelope. On A and P, λa=α^-nλ exactly, because D∞⊂Ha. Choose a≥8n^4; then κ+2/(a+1)<αβ. Every actual target item meeting A has radius <3, while all source items meeting P have radius >4. The same actual target j therefore receives only earlier source items, and the center-displacement guard makes each source bQi capture yj. Its one original eligibility bounds rowμa≤τ, including the radial and Y restrictions.

Writing ca=α^-n≥1, ambient energy scales by ca², degrees by ca and target mass by ca. Hence the finite-a deficit is at least

\[
\frac{c_a^2E_{\rm good}-\tau c_a a_*}{c_aHM_0}
\ge\frac{E_{\rm good}-\tau a_*}{HM_0}
\ge2^{-31}nd.
\]

This holds for every a≥8n^4 and every original finite history at all sufficiently large n. No assumption on how μ is distributed, and no convergence of selectors or histories, is used.

## 8. The new scalar scale makes τ negligible and yields the new order

Set δ=2^-31nd and τ=√[A0(B0−δ)]. The original cap integral estimates, with x∈[L/2,L], give B0>8n eventually, so δ<B0. Also

\[
A_0\le e^x=ne^{-s},\qquad B_0\le nL,
\quad\tau\le n\sqrt L\,e^{-s/2}.
\]

Therefore

\[
\frac\tau{M_0}
\le100L^{7/6}e^{-s/4}\longrightarrow0.
\]

Thus the condition τ/M0≤cE/4 in Section 7 holds eventually. Insert the resulting deficit into the accepted original finite cover/load formula and let a→∞ only in its uniform scalar bounds. This yields (1), preserving general L1 inputs and every positive common radius.

Finally x=L(1−d²), so

\[
B_0\le1+n\log t
\le T+2nd^2+1
\]

eventually, using −log(1−d²)≤2d². Further,

\[
A_0+2\sqrt{A_0B_0}
\le ne^{-s}+2n\sqrt L e^{-s/2}=o(nd).
\]

Consequently F(A0,B0−δ)≤T−2^-31nd+o(nd)≤T−2^-32nd eventually, proving (2).

## 9. Scope and quantitative status

The proposed improvement uses a different x and a second-moment column cutoff. It does not use the R035 global bad-source/M0 argument, which would fail at this scale. The fixed-face coordinate, target interval endpoints, monotonicity in target r, and the column-deletion energy direction are all explicit above.

The theorem claims only a common sufficiently-large-n threshold. No finite threshold, global parameter optimum, or claim that C2=144 is sharp is made. The independent one-dimensional moment audit may use a smaller constant, but that is not needed for the rational constants in this proof. No finite numerical batch can replace the full operator bridge in Section 7.

