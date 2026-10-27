# R051 draft · A flexible exact budget family for the bounded inverse-column mechanism

Status: generic full-operator proof submitted for independent review. It does not certify a numerical search before that review. All parameters below are fixed independently of n; every asserted eventual threshold may depend on the complete fixed tuple. Old sources remain unchanged. No n=7 work, production search, MATLAB, or Git operation.

## 1. Complete parameter domain and rational budget

Fix real parameters

\[
a>3,\qquad s>2a,\qquad
0<d\le1/16,\qquad d<1/(32a),\qquad
\sigma=sd^2<1/2.
\]

The letter a here is the **soft-weight exponent**, not the original small-cube dilation; the latter is written mathfrak a below. Fix errors

\[
0<\varepsilon_J,\varepsilon_A,\varepsilon_D<1,
\qquad\varepsilon_S,\varepsilon_T>0.
\]

These upper bounds on εJ,εA,εD are part of the domain, not optional: otherwise denominators or signs can produce spurious positive budgets.

Set

\[
P_0=\frac{3579}{4096},\quad q_{\min}=1-3d/4,
\]

\[
u_{\rm lo}=\frac{P_0}{2(1+\varepsilon_J)},\qquad
u_{\rm hi}=\frac1{P_0(1-d/2)(1-\varepsilon_J)},
\]

\[
M_* =\max\left(4,\frac{2u_{\rm hi}}{q_{\min}}\right),
\qquad
b_*=(1-3/a)u_{\rm lo}-\varepsilon_S-\varepsilon_T.
\]

Require b*>0. Define the exact mass, deficit and net-cost coefficients

\[
m_*=(1-\varepsilon_A)\frac{d}{2(1+3d/4)},\qquad
\Gamma=\frac{m_*b_*}{M_*},\qquad
\gamma=(1-\varepsilon_D)\Gamma,
\]

\[
c_{\rm raw}=\gamma-\frac{\sigma}{1-\sigma}. \tag{1}
\]

For rational input parameters all the displayed budget quantities are exact rational numbers; no exponential approximation is needed to check them.

**Generic theorem.** With L=logn, T=nL−nlogL and

\[
x=(1-\sigma)L,\quad t=n/x,\quad\beta=1+1/t,
\quad A_0=\beta^n,\quad B_0=B_n(t),
\]

every tuple in the domain satisfies, for all sufficiently large n,

\[
C_n\le F(A_0,B_0-\gamma n),\qquad F(A,B)=(\sqrt A+\sqrt B)^2. \tag{2}
\]

If c_raw>0, then every fixed 0<c<c_raw gives

\[
C_n<T-cn \tag{3}
\]

for all sufficiently large n. This is an original centered-cube maximal-operator upper bound for general L1 inputs and all positive common radii, not merely a model objective. A numerical tuple can invoke it only after the entire domain and its strict rational inequalities are checked.

## 2. Geometry and all eventual guards

Let ε0=2^-20, h=ε0/(Ld), ζ=n^-3, κ=β−ζ, θ=exp(−aLd²), and k=nhθ. Targets have r=1+u with u∈I=[d/4,3d/4]. Sources have ρ=tq, q=1−v with v∈I. Matching is |r−q−d|<h.

Use the same ambient corner measure

\[
d\lambda(z)=2^{-n}m(z)^{-n}\mathbf1_{\mathcal D}(z)dz,
\quad m(z)=\max(1,\|z\|_\infty),
\]

\[
\mathcal D=\{\|z\|_\infty<t:\exists i,z_i<1-\beta m(z)\}.
\]

Before cap restriction, a face at radius R>1 has measure dR/(2R) and n−1 independent Uniform[−R,R] free coordinates.

Since 0<σ<1/2, x lies between L/2 and L, x→∞ and t∼n/L→∞. Consequently all of the following hold eventually for each fixed tuple:

- κ>1 because 1/t≫n^-3; β<2.
- h<d/8, ζt<h/2 and Ld≥1.
- Every target r<2 and every source ρ>4r>4.
- q≥qmin>1/2 and the baseline q0=r−d+h remains in a fixed interval near 1, even if it lies h beyond the source-band endpoint.
- With R(q)=q+d−h and g=d−h+ζtq, one has R>1, 0<g≤d≤1/16 and ρ≥4R.

No threshold uniform over an unfixed continuous family of tuples is asserted.

## 3. Typical target product and the sharper individual-coordinate bounds

For the exact inverse row set

\[
\ell_R(z)=\min(1,(R+\kappa\rho-|z|)/(2R)),\quad
J(w;q)=\frac1{2\rho}\int_{[-\rho,\rho]\cap(w-\kappa\rho,w+\kappa\rho)}\ell_R(z)^{-1}dz.
\]

Its closed form is J=1+[2B−A(T)]/(2ρ), where T=(|w|−(κ−1)ρ)+,

\[
A(T)=\begin{cases}2R\log(1+T/(2R-g)),&T\le g,\\
2R\log(2R/(2R-g))+T-g,&T\ge g,\end{cases}
\quad B=A(g)-g.
\]

The local proof uses only the guards of Section 2, not the particular soft exponent or scalar buffer. At q0=r−d+h, R=r. For a free target coordinate W uniform on [−r,r], Fubini gives EJ=1 exactly. The elementary conditional moment bounds are

\[
-6g^3/\rho^2\le E\log J\le0,\qquad
\operatorname{Var}(\log J)\le12g^3/\rho^2,
\]

while the pinned coordinate ±r contributes O(g/ρ). Thus the n-coordinate logarithmic sum tends to zero in probability on every target face, uniformly in r, since nd³/t²=O_d(L²/n)→0.

Let A be the target annulus intersected with D and the two conditions

\[
|\sum_i\log J(w_i;q_0)|\le1/8,
\qquad K_w:=\#\{|w_i|>r-2d\}\le4nd.
\]

The second condition also has probability 1−o(1): on a fixed face Kw=1+Bin(n−1,2d/r), with mean≤1+2nd and variance≤2nd. The target corner-captured fraction is exponentially small in n. Therefore for the arbitrary fixed εA>0 of Section 1,

\[
\lambda(A)\ge(1-\varepsilon_A)n
\log\frac{1+3d/4}{1+d/4}
\ge m_*n \tag{4}
\]

eventually.

For the full matching interval, the locally Lipschitz function logJ obeys, a.e. in q,

\[
|\partial_q\log J(w;q)|
\le\frac{64}{\rho}[d^2+\mathbf1_{\{|w|>r-2d\}}].
\]

This follows by differentiating the displayed A/B form; capture onset is a kink, so the argument uses absolute continuity rather than everywhere classical differentiability. Summing over coordinates and integrating over |q−q0|≤2h gives at most 1280Ldh=5/4096. Hence on A, uniformly over matched source radii,

\[
|\sum_i\log J(w_i;q)|\le\frac{517}{4096}.
\]

The inequalities e^-z≥1−z and e^z≤1/(1−z), valid here with z=517/4096<1, prove the exact rational enclosure

\[
P_0\le\prod_iJ(w_i;q)\le1/P_0. \tag{5}
\]

In addition, the closed form gives |J−1|≤C d/ρ uniformly for all target-coordinate values |w|≤r and all matched or baseline q. Since ρ→∞, the arbitrary fixed εJ gives eventually the pointwise bound

\[
1-\varepsilon_J\le J(w;q)\le1+\varepsilon_J. \tag{6}
\]

This does not require a union bound over n coordinates. A separate coarse geometric J≥1/2 remains available for the inverse-probability estimate below, even if εJ itself is larger than 1/2.

## 4. Inverse rows, soft retention and exact degree budgets

Define the source product and weight

\[
G(z)=\prod_i\ell_R(z_i)>0,\qquad
b(z)=\min(1,\theta/G(z))\in[0,1].
\]

The full inverse-weighted source shell coefficient at target w is exactly

\[
U(w,\rho)=\frac{\prod_jJ(w_j;q)}{2\ell_*}
\sum_i\frac{N_i}{J(w_i;q)},
\quad\ell_*=1-g/(2R),
\quad N_i=1+\mathbf1_{\{|w_i|<(\kappa-1)\rho\}}.
\]

Since 1≤Ni≤2, 1−d/2≤ell*≤1, and (5)–(6) hold,

\[
u_{\rm lo}n\le U(w,\rho)\le u_{\rm hi}n. \tag{7}
\]

For each captured source face, normalize the inverse capture density proportional to 1_capture/G. Its free coordinates are independent with densities proportional to 1_capture/ell_R and normalizers J≥1/2. Using ell_R≥3/4,

\[
E_{\rm inv}[-\log\ell_R(Z_i)]
\le2E_{\rm full}[(-\log\ell_R)/\ell_R]
\le8g^2/(9R\rho)\le d^2/\rho.
\]

The pinned coordinate contributes at most d. Thus on each face and their inverse-row mixture,

\[
E_{\rm inv}[-\log G]\le d+nd^2/\rho
\le d+2Ld^2\le3Ld^2.
\]

Markov in this inverse-row probability gives P_inv(G<θ)≤3/a. Since a>3, the full soft shell coefficient is between

\[
(1-3/a)u_{\rm lo}n\theta
\quad\text{and}\quad u_{\rm hi}n\theta. \tag{8}
\]

Every target matching interval has ordinary v-width at least h and at most 2h. Its source log-width is therefore between h and 2h/qmin. The raw soft row is bounded below by (1−3/a)u_lo k and above by 2u_hi k/qmin.

Because b≤1, removing the original corner-captured source part costs at most its unweighted full-shell mass. That mass is O(nd e^-Ld/32), using n≥2 and x≥L/2. The ratio to k is

\[
O\left(\frac dh\exp[-(d/32-ad^2)L]\right)\longrightarrow0
\]

by the strict condition d<1/(32a). Thus for the arbitrary fixed εS>0, the final row lower bound is

\[
R_\lambda(w)\ge[(1-3/a)u_{\rm lo}-\varepsilon_S]k. \tag{9}
\]

For every source, the actual target column is bounded by the full matched column

\[
C_{\rm raw}(z)\le(4nh/R)G(z).
\]

Indeed target side lengths decrease with target radius, r≥R, their minimum is at least 1/2, and matched log-width≤2h/R. Therefore

\[
b(z)C_{\rm raw}(z)\le4k.
\]

Together these give the exact maximum-degree budget

\[
\max(R_\lambda,C_\lambda)\le M_* k. \tag{10}
\]

No source concentration or product-moment cutoff is used. The potentially large source variation is absorbed by the bounded position-dependent weight and the recomputed row scale.

## 5. Complete choice of τ and the same-history energy deficit

Set δ=γn and τ=√[A0(B0−δ)]. The original cap estimates with L/2<x<L give B0>8n eventually. Also Γ<d/16≤1/256 from the displayed positive budgets, so δ<B0. Further

\[
A_0\le e^x=ne^{-\sigma L},\quad B_0\le nL,
\quad\tau\le n\sqrt L e^{-\sigma L/2}.
\]

Consequently

\[
\frac\tau k\le\frac{\sqrt L}{h}
\exp[-(s/2-a)d^2L]\longrightarrow0
\]

because s>2a. Thus τ≤εT k eventually. By (4), (9)–(10), the actual energy deficit lemma would give

\[
\Delta\ge\frac{\lambda(A)[(1-3/a)u_{\rm lo}-\varepsilon_S-\varepsilon_T]k}{M_*k}
\ge\Gamma n>\delta. \tag{11}
\]

It remains to check that the displayed energy inequality really applies to every original finite history, not merely the ambient model; the next section performs that check.

## 6. Every original μ, general L1 witnesses and every positive scale

Use mathfrak a≥8n^4 for the old small-cube dilation and put

\[
\alpha=\frac{\mathfrak a}{\mathfrak a+1},\qquad
c_{\rm cap}=\frac{\mathfrak a\beta-1}{\mathfrak a+1}.
\]

At any load support point, normalize the last minimum positive-contribution radius to 1 and reflect the last center into the positive orthant. The actual uncaptured increment measure μa is dominated by the old corner envelope λa. On the new target and source sets, λa=α^-nλ exactly, since

\[
\alpha-c_{\rm cap}m-(1-\beta m)
=\frac{(\beta+1)m-1}{\mathfrak a+1}>0.
\]

Every actual target meeting A has radius<3, from its uncaptured-center inequality and ||w||<2. Every source item has radius>4 because its increment coordinate has ||z||=ρ>4. Sources therefore precede the same actual target j. The guard

\[
\kappa+2/(\mathfrak a+1)<\alpha\beta
\]

holds for every mathfrak a≥8n^4, so the true centers are captured. Multiplying included source increments by b≤1 preserves j's eligibility budget τ. Thus weighted actual rowμa≤τ on μa-a.e. target, without any distribution assumption on μa.

The bipartite energy identity uses the single missing measure λa−μa and disjoint target/source sets. Ambient energy scales by α^-2n, row and column degrees by α^-n, and target mass by α^-n. Because α^-n≥1, this scaling can only improve the lower bound in (11). Hence every original finite history has total cap deficit at least Γn, and in particular at least δ=γn, uniformly for all mathfrak a≥8n^4 once n is large enough for this fixed tuple.

The remaining captured budget is τ and the terminal self-weight at most mathfrak a^-n. The accepted original cover and witness argument therefore gives

\[
C_n\le\left(1+\frac{(\beta+1/\mathfrak a)^n}{\tau}\right)
\left[\tau+\mathfrak a^{-n}+\alpha^{-n}J_{\mathfrak a}(t)-\delta\right].
\]

Let mathfrak a→∞ only in this family of uniform scalar inequalities. This proves (2). Original finite covers, actual increments, general L1 input witnesses and the supremum over all positive common radii are unchanged; weights are only a post-selection energy tool. No selector convergence or independently chosen histories are assumed.

## 7. Full scalar budget and exact old/new comparison

Since x=(1−σ)L,

\[
B_0\le1+n\log t
=T-n\log(1-\sigma)+1
\le T+\frac{\sigma}{1-\sigma}n+1.
\]

Also A0+2√(A0B0)≤ne^-σL+2n√L e^-σL/2=o(n), as σ>0 is fixed. Thus the full new root cost, not just its cap component, satisfies

\[
F(A_0,B_0-\gamma n)
\le T-c_{\rm raw}n+o(n).
\]

Every fixed c<c_raw has positive slack, proving (3) when c_raw>0. The endpoint c=c_raw is not automatically supplied by this asymptotic upper envelope.

At the same t, the exact old/new difference is

\[
F(A_0,B_0)-F(A_0,B_0-\gamma n)
=\gamma n+2\sqrt{A_0}
[\sqrt{B_0}-\sqrt{B_0-\gamma n}]>\gamma n.
\]

Against other accepted t choices or explicit thresholds, use the minimum of all bounds whose hypotheses are met. Rational budget improvements among tuples do not assert common finite-dimensional validity, a uniform threshold over the continuous parameter domain, or parameter global optimality.

## 8. Certification interface and what is not proved

For each proposed fixed tuple, an exact verifier must check the full domain in Section 1, including b*>0, and compute gamma and c_raw with outward or exact arithmetic. If c_raw>0, choose a fixed exact c strictly inside (0,c_raw). After the generic theorem is independently accepted, these checks certify an eventual original-operator bound. Before that acceptance they remain conditional candidates.

The εJ,εA,εD upper bounds prevent false sign reversals; positivity of b* prevents a negative row budget being disguised by another negative factor. If any domain condition fails, this theorem supplies no claim even if a reduced spreadsheet expression happens to be positive.

The proof has only a common sufficiently-large threshold for each fixed tuple. A tuple chosen as a function of n, a whole adaptive continuous search, or a finite explicit n threshold requires additional uniform analysis. No such statement is implicit in the budget family.

