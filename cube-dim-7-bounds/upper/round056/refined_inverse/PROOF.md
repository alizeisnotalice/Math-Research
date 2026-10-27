# R056 draft · Refined inverse rows, exact column factors, and a finite general bound

Status: full proposed operator proof for independent review. The intended claim is not accepted merely because its rational budget is positive. This source includes the actual-history, general-L1/all-scale bridge and an explicit finite guard system. It does not modify any earlier accepted source. No n=7 work, large-CPU search, shared MATLAB, or Git operation.

## 1. Fixed data and proposed theorem

Let L=logn and require

\[
L\ge L_0=2^{20},\quad d=3/10,\quad
\sigma=7/100,\quad x=(1-\sigma)L=93L/100,\quad t=n/x,
\quad\beta=1+1/t.
\]

Fix

\[
\chi=7/200,\quad\epsilon_0=2^{-24},\quad
h=\frac{\epsilon_0}{Ld},\quad\zeta=n^{-3},\quad
\kappa=\beta-\zeta,\quad\theta=e^{-\chi x},\quad k=nh\theta.
\]

The source band is q=ρ/t∈[qmin,qmax]=[31/40,37/40]. Targets have radius

\[
r\in[r_{\rm lo},r_{\rm hi}],\qquad
r_{\rm lo}=43/40+2^{-20},\quad r_{\rm hi}=49/40-2^{-20}.
\]

Set all εJ,εN,εRet,εS,εT,εD equal to e=1/1024, and set εA=1/64. Define the exact rational quantities

\[
P_0=65275/65536,\quad
U_{\rm lo}=P_0(1-e)/(1+e),\quad
U_{\rm hi}=(1+e)/[P_0(1-e)],\quad M=21/8.
\]

Let zA=(rhi−rlo)/rlo and use the four-term lower log bound

\[
m_A=(1-\epsilon_A)
[z_A-z_A^2/2+z_A^3/3-z_A^4/4],
\]

\[
b_*=(1-e)\frac{2U_{\rm lo}}{q_{\max}}-2e,
\qquad\gamma=(1-e)\frac{m_Ab_*}{M}.
\]

**Proposed theorem.** For every integer n with logn≥2^20,

\[
C_n\le F(A_0,B_0-\gamma n)
<n\log n-n\log\log n-\frac n{40}, \tag{1}
\]

where A0=β^n, B0=B_n(t), F(A,B)=(√A+√B)^2. The exact budget checker gives γ−7/93>1/40, but the complete theorem requires all proofs below. There is no search-optimality or public-priority claim.

## 2. Ambient geometry, the target baseline, and local d=0.3 guards

Use the original corner envelope

\[
d\lambda=2^{-n}m(z)^{-n}\mathbf1_{\mathcal D}(z)dz,
\quad m(z)=\max(1,\|z\|_\infty),
\]

\[
\mathcal D=\{\|z\|_\infty<t:\exists i,z_i<1-\beta m(z)\}.
\]

Before restriction to D, each face at radius R>1 has measure dR/(2R); its other n−1 coordinates are independent Uniform[−R,R]. Full shell mass is n dR/R.

For a target radius r define the matching center qc=r−d and the two endpoints q±=qc±h. Matching means |q−qc|<h. The **log-product baseline is qb=q+**, not qc. At qb,

\[
R(q)=q+d-h=r,
\]

so the exact baseline identity EJ=1 is available. Since h≤2^-20, the trimmed target interval ensures the entire matching interval [q−,q+] lies inside the source band. Its ordinary width is exactly 2h; endpoint losses from the untrimmed family do not occur.

For any matching q put

\[
\rho=tq,\quad R=q+d-h,\quad
a_0(q)=(\kappa-1)\rho=(1-\zeta t)q,
\quad g=R-a_0=d-h+\zeta\rho.
\]

The finite guards below imply q≥1/2, R≥10749/10000>1, ρ≥4r, 0<g≤d≤min(1/3,R/2), h≤d/8 and κ>1. In particular the d≤1/16 assumption of older budget statements is **not** blindly reused.

Define

\[
\ell_R(z)=\min(1,(R+\kappa\rho-|z|)/(2R)),
\quad J(w;q)=\frac1{2\rho}\int_{[-\rho,\rho]\cap(w-\kappa\rho,w+\kappa\rho)}\ell_R(z)^{-1}dz.
\]

The independently checked local lemma, valid for d≤1/3 under the displayed guards, has the following statements. Its derivation uses

\[
F(R,s)=-2R\log(1-s/(2R))-s,
\quad m=(|w|-a_0)_+,\quad s_0=(g-m)_+,
\]

\[
J=1+[F(R,g)+F(R,s_0)-m]/(2\rho).
\]

The bounds Fs≤d, |FR|≤d² and F≤(3/10)d² give, a.e. in q,

\[
|\partial_q\log J|\le(64/\rho)[d^2+\mathbf1_{\{|w|>r-2d\}}]. \tag{2}
\]

The capture-onset kink is handled by local Lipschitz continuity and absolute continuity, not a claim of everywhere differentiability. At R=r the exact mean EJ=1 and the conditional-face bounds are

\[
-6g^3/\rho^2\le E\log J\le0,
\qquad\operatorname{Var}(\log J)\le12g^3/\rho^2.
\]

The pinned coordinate contributes a deterministic term in [−g/ρ,0]. These constants require R≥1,g≤R/2,ρ≥4R, not d≤1/16. They are re-established for this expanded domain in the independent local lemma at `audit/round056/one_coordinate/LEMMA.md`.

## 3. Target set: log product and normalized source-face count

For fixed w with r=||w|| define

\[
N_i(q)=1+\mathbf1_{\{|w_i|<a_0(q)\}},\quad
T_N(q)=\sum_iN_i(q),\quad
\ell_*(q)=1-g/(2R)=\frac{R+a_0}{2R}.
\]

On a fixed maximal target face the pinned coordinate has Ni=1, since a0(q)<r. The remaining n−1 coordinate indicators are independent Bernoulli with parameter a0(q)/r. Thus

\[
E_\pm=n+(n-1)a_0(q_\pm)/r,
\qquad\operatorname{Var}T_N(q_\pm)\le n/4.
\]

Let A be the target annulus intersected with D and these three conditions:

\[
|\sum_i\log J(w_i;q_b)|\le1/256,
\]

\[
|T_N(q_-)-E_-|\le e n/4,
\quad |T_N(q_+)-E_+|\le e n/4.
\]

For d=0.3, K_w=#{|wi|>r−2d}≤n≤4nd is automatic. Integrating (2) from qb over length≤2h yields

\[
|\sum_i\log J(w_i;q)-\sum_i\log J(w_i;q_b)|
\le1280Ldh=5/65536.
\]

Hence the full matching product obeys

\[
P_0\le\prod_iJ(w_i;q)\le1/P_0,\quad
P_0=1-(1/256+5/65536)=65275/65536. \tag{3}
\]

### Uniform count calibration, not independent face worst cases

Both TN(q) and ell*(q) increase with q: a0 is increasing, and ell*'=(1−ζt)(d−h)/(2R²)>0. Since R(q−)=r−2h and R(q+)=r, the two endpoint concentration events imply, for every matching q,

\[
\frac{T_N(q)}{2n\ell_*(q)}
\ge1-[2h+1/n+e/4].
\]

For the upper bound, use

\[
\frac{a_+}{r}-\frac{a_-}{r-2h}
=\frac{2h(1-\zeta t)(d-h)}{r(r-2h)}\le2hd,
\]

which gives TN(q)/(2nell*)≤1+2hd+e/4. Therefore the finite gate 3h+1/n≤e/2 ensures

\[
1-e\le T_N(q)/(2n\ell_*(q))\le1+e. \tag{4}
\]

The same two events control all q; no independent maximization at different positions is summed.

### Quantified good mass

Let N(L)=L^8/8!≤n. Uniform baseline bounds give mean absolute value at most

\[
M_L=\frac{6d^3L^2}{q_{\min}^2N(L)}+
\frac{dL}{q_{\min}N(L)}
\]

and variance at most

\[
V_L=\frac{12d^3L^2}{q_{\min}^2N(L)}.
\]

If ML≤1/512 and VL≤2^-26, the baseline bad probability is≤1/256=εA/4. Each count bad probability is≤4/(e²n)≤4/(e²N(L)), also≤εA/4 by the gates. On any single target face the original corner-captured fraction is≤(3/4)^(n−1): a negative maximal face contributes zero, and a positive maximal face contributes ell_t(r)^(n−1), with t≥4 and r<49/40 ensuring ell_t(r)<3/4. The gate n≥75 makes this per-face fraction≤εA/4. Union of these four bad events therefore has probability≤εA on every full target face and radius. Consequently

\[
\lambda(A)\ge(1-\epsilon_A)n\log(r_{\rm hi}/r_{\rm lo})\ge m_A n. \tag{5}
\]

The log lower bound in mA is the four-term alternating Taylor lower bound at 0<zA<1.

## 4. Refined full inverse row and exact column factor

Define the full source band P0={z:rho=||z||infty, rho/t∈[qmin,qmax]}. On all of P0, including points outside the corner cap, define

\[
G(z)=\prod_i\ell_{R(\rho/t)}(z_i),\qquad
b(z)=\min(1,\theta/G(z)).
\]

Set b=0 off P0. This extension is made before any Dcap restriction; the full-row inverse laws below use the whole band and are not cap-conditioned. Let P=D∩P0 and let A be the target set from Section 3. For r=||w||infty and rho=||z||infty the final kernel is exactly

\[
K(w,z)=\mathbf1_A(w)\mathbf1_P(z)
\mathbf1_{\{|\rho/t-(r-d)|<h\}}
\mathbf1_{\{\|z-w\|_\infty<\kappa\rho\}}\,b(z).
\]

Its ambient degrees are Rlambda(w)=∫K(w,z)dλ(z) and Clambda(z)=∫K(w,z)dλ(w). The raw soft row is the corresponding integral over the whole P0 against 2^-n rho^-n dz, before deleting P0\D in Section 6. A and P are disjoint: r≤rhi<2 whereas the finite source-radius gate gives rho≥4rhi>4.

The same finite guards ensure the pointwise bound J∈[1−e,1+e] for every coordinate: the local formula gives |J−1|≤2d/ρ, and the explicit gate 2dL/[qmin N(L)]≤e suffices. The inverse shell coefficient is

\[
U(w,\rho)=\frac{\prod_jJ(w_j;q)}{2\ell_*}
\sum_i\frac{N_i(q)}{J(w_i;q)}.
\]

Combining (3), (4), and individual J bounds yields

\[
U_{\rm lo}n\le U(w,\rho)\le U_{\rm hi}n. \tag{6}
\]

The pin factor is canceled by the normalized count, rather than separately losing one half in the lower estimate.

For a fixed source z, put ai=κρ−|zi|>0. At target radius r, normalized side lengths are ell_i=min(1,(1+ai/r)/2). If ai≥r then ell_i=1 and the target face count Ni≤2. If ai<r, exactly one target face is captured and ell_i>1/2. Thus always Ni/(2ell_i)≤1. The exact full target-shell coefficient is therefore at most n times its n-side product, rather than 2n times that product. Since r≥R and each side length decreases with r,

\[
C_{\rm raw}(z)\le(2nh/R)G(z),
\quad G(z)=\prod_i\ell_R(z_i).
\]

With b(z)=min(1,θ/G(z)), the final weighted column satisfies

\[
C_\lambda(z)\le2k. \tag{7}
\]

## 5. Per-source-face inverse Chebyshev, then face mixing

Let φ=−logell_R, Qpin=−logell*. The exact full-source free-coordinate identities are

\[
E_{\rm full}[\phi/\ell]=R Q_{\rm pin}^2/\rho,
\quad E_{\rm full}[\phi^2/\ell]=2R Q_{\rm pin}^3/(3\rho).
\]

On a fixed captured source face the pin contributes Qpin deterministically; its other n−1 coordinates are independent under their inverse-capture laws, with normalizers J≥1−e. Hence

\[
E(V\mid\text{face})\le Q_{\rm pin}
+x\frac{R}{q(1-e)}Q_{\rm pin}^2,
\]

\[
\operatorname{Var}(V\mid\text{face})
\le x\frac{2R}{3q(1-e)}Q_{\rm pin}^3,
\qquad V=-\log G.
\]

The mean generally depends on the face. We do not assert that the mixture variance equals the sum of within-face variances.

For h≤1/10000, R≥10749/10000, g≤d, and R/q≤43/31. Thus

\[
Q_{\rm pin}\le\log(3583/3083)
\le Q:=u-u^2/2+u^3/3,\quad u=500/3083.
\]

The exact rational checks give Q<1/6 and

\[
\frac{43}{31(1-e)}Q^2<4/125,
\qquad\frac{2}{3}\frac{43}{31(1-e)}Q^3<1/300.
\]

Therefore every captured face has mean≤Q+4x/125 and variance≤x/300. For x≥1000, Q/x≤1/2000, so the threshold χx exceeds this mean by at least

\[
x[\chi-4/125-1/2000]=x/400.
\]

Chebyshev on **each face separately** gives

\[
P(V>\chi x\mid\text{face})\le1600/(3x)\le e.
\]

Mixing these probability inequalities preserves the bound e, without a claim about mixture variance. This is done in the full inverse-capture source law before corner-cap deletion. Thus soft retention gives shell coefficient≥(1−e)Ulo nθ, while its upper coefficient is≤Uhi nθ.

Since the matching interval has exact width 2h inside [qmin,qmax], its log-width is between 2h/qmax and 2h/qmin. Before cap deletion, rows are therefore bounded below by (1−e)2Ulo k/qmax and above by 2Uhi k/qmin.

## 6. Original cap, τ, and the new exact degree budget

The original source corner-captured coefficient is (n/2)ell_t(ρ)^(n−1). Since q≤37/40 and n≥75,

\[
\frac{n-1}{2}\left(\frac1\rho-\frac1t\right)
= (1-1/n)x\frac{1-q}{2q}
\ge x/25.
\]

The difference 1/25−χ equals 1/200. Since b≤1, the removed weighted row is bounded by the full unweighted cap mass, and hence by

\[
nd e^{-x/25}\le k\,(d/h)e^{-x/200}\le ek. \tag{8}
\]

The finite gate for the last inequality is recorded below. The old rough exp(−Ld/32) estimate would be too weak here and is not used.

Choose τ=√[A0(B0−γn)]. The accepted cap estimates give B0>8n, while γ<1; so the budget is positive. A0≤ne^-σL and B0≤nL imply

\[
\frac\tau k\le\frac{\sqrt L}{h}
e^{-49L/20000}\le e. \tag{9}
\]

All final weighted rows are thus at least [(1−e)2Ulo/qmax−e]k; after paying τ their excess is at least b*k. Exact rational checks give

\[
2U_{\rm hi}/q_{\min}<21/8,
\]

and columns are≤2k by (7), so maximum ambient degree is≤Mk. The same-history energy deficit is at least mA b* n/M, which strictly exceeds γn because εD=e>0.

## 7. Finite guard table, valid for all L≥L0

Let N(L)=L^8/8! and use x=93L/100. It is enough to check the following rational inequalities at L=L0:

1. N≥75, x≥1000, N/L≥4, qmin N/L≥4rhi.
2. h≤min(2^-20,d/8,1/10000), and 1/[xN²]≤h/2. These give band containment, g≤d and κ>1.
3. 2dL/(qmin N)≤e; ML≤1/512; VL≤2^-26.
4. 4/(e²N)≤εA/4; 3h+1/N≤e/2; (3/4)^74≤εA/4.
5. 1600/(3x)≤e and Q/x≤1/2000.
6. The following Taylor-eight bounds are at most e:

\[
\frac{8!d^2}{\epsilon_0(93/20000)^8L^7},
\qquad
\frac{8!d}{\epsilon_0(49/20000)^8L^{13/2}}.
\]

They bound (d/h)e^-x/200 and (√L/h)e^-49L/20000, respectively, using exp(y)≥y^8/8!. At L0, √L0=2^10 makes every expression rational.

All gates extend to larger L: their nonconstant upper bounds are negative powers of L, or positive sums of such powers. Lower bounds N/L, x and N are increasing. The ratio [1/(xN²)]/h is a negative power L^-16. Thus no finite sampling is used for L>L0.

The scalar final envelope is controlled by the rational value at L0 of

\[
P(L)=\frac{8!}{(\sigma L)^8}
+\frac{2\sqrt L\,8!}{(\chi L)^8}
+\frac{8!}{L^8},
\]

whose three terms are decreasing powers of L. Exact arithmetic verifies

\[
P(L_0)<\gamma-7/93-1/40.
\]

The source `verify.py` independently reconstructs all these rational checks, the log Taylor bounds, mass coefficient, degree coefficient and net margin. It is one fixed theorem tuple, not a production search.

## 8. All original histories, general L1 inputs and all positive scales

Use the original dilation mathfrak a≥8n^4, α=mathfrak a/(mathfrak a+1), and ccap=(mathfrak aβ−1)/(mathfrak a+1). At any load support point normalize the last minimum positive-contribution radius to 1 and reflect its center into the positive orthant. The actual uncaptured-increment measure μa obeys μa≤λa, and λa=α^-nλ exactly on the new target and source sets, since D∞⊂Ha.

All actual target items meeting A have radius<3 from their uncaptured-center inequality and ||w||<2. For an actual source increment point z=y/(mathfrak a+1), the common support 0∈mathfrak a Qi and y∈Qi give ||z||infty<ri, so ri≥m(z). The source-radius gate gives rho≥4rhi>4, hence every actual source item has radius>4. Thus every source precedes the same actual target j. The guard κ+2/(mathfrak a+1)<αβ makes the actual centers captured. Since b≤1, this same j's eligibility pays the weighted row≤τ, with no inverse factor paid by actual μ.

The bipartite identity uses the single deficit λa−μa and disjoint A/P sets. Energy scales by α^-2n, degrees and target mass by α^-n, and α^-n≥1 can only improve the deficit bound. Hence every original finite history has missing cap mass at least mA b*n/M>γn, uniformly for mathfrak a≥8n^4. The captured budget remains τ and the final self weight at most mathfrak a^-n.

The original finite cover and witness argument therefore yields

\[
C_n\le\left(1+\frac{(\beta+1/\mathfrak a)^n}{\tau}\right)
\left[\tau+\mathfrak a^{-n}+\alpha^{-n}J_{\mathfrak a}(t)-\gamma n\right].
\]

Let mathfrak a→∞ only in these uniform scalar inequalities. This proves the first bound in (1), without requiring convergence of histories or selectors. It retains the original general-L1 witnesses and supremum over every positive common radius; no truncated operator or independent witness family is substituted.

For the complete scalar cost, B0≤T+(7/93)n+1 and A0≤ne^-σL. Thus

\[
\frac{F(A_0,B_0-\gamma n)-T}{n}
\le-\gamma+7/93+e^{-\sigma L}
+2\sqrt L e^{-\chi L}+e^{-L}
\le-\gamma+7/93+P(L)<-1/40.
\]

This proves the strict final claim. It is not inferred from a cap-component improvement alone.

## 9. Evidence and limitations

`verify.py` is a small standard-library Fraction verifier for the exact fixed tuple and finite guard bases. Analytic monotonicity extends those bases to every L≥L0. The source code does not optimize parameters, replay any old grid, or certify a global optimum. Independent numerical kernel diagnostics are supporting evidence only, not proof of untested domains.

The actual one-coordinate lemma explicitly allows d≤1/3 with R≥1 and g≤R/2; the proof does not extend an old d≤1/16 theorem by name alone. Face-count concentration is under the full target-face law; inverse Chebyshev is under each captured full source-face law; cap restriction is paid afterwards. These distinct probability spaces must not be interchanged.

