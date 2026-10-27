# R061 draft · Disjoint-layer energy payment on one actual history

Status: proposed full operator result for independent review. The exact R062 record is a conditional scalar certificate until the geometric and actual-history theorem below is accepted. This proof uses one common normalization and one global selection threshold; it never makes separate histories for different layers. No n=7 work, production grid replay, MATLAB, or Git operation.

## 1. Direct-sum lemma for a single deficit measure

Let λ be finite, μ≤λ, and D0=λ−μ≥0. Suppose Borel pairs (Aj,Pj), j=1,...,J, have pairwise disjoint unions Bj=Aj∪Pj, up to λ-null sets; in particular Aj∩Pj is null. Let Kj∈[0,1] be a Borel kernel supported in Aj×Pj, with

\[
\int K_j(w,z)d\mu(z)\le\tau
\quad\text{for }\mu\text{-a.e. }w\in A_j.
\]

The same τ is used for every layer. Suppose its ambient row lower degree is rj kj, and both row and column upper degrees are at most Mj kj, with kj>0. Then

\[
D_0(B_j)\ge
\lambda(A_j)\frac{r_j-\tau/k_j}{M_j},
\qquad
D_0(H)\ge\sum_j
\lambda(A_j)\frac{r_j-\tau/k_j}{M_j}. \tag{1}
\]

Here H contains all supports; negative local bounds may be replaced by zero.

**Proof.** For the same actual μ and λ,

\[
E_{\lambda\lambda,j}-E_{\mu\mu,j}
=E_{D_0|A_j,\lambda|P_j}+E_{\mu|A_j,D_0|P_j}
\le M_jk_j[D_0(A_j)+D_0(P_j)].
\]

Actual row payment gives Eμμ,j≤τμ(Aj)≤τλ(Aj), while ambient energy is at least λ(Aj)rjkj. This proves the local inequality. Sum D0(Bj), using actual disjointness of these point supports, to obtain (1). ∎

The lemma does not partition the original item labels. If one actual increment crosses several source bands, its measure is split pointwise; the sum of its restricted masses is at most its one original weight. Nothing is duplicated. The same is true if one actual target increment meets multiple Aj. All local inequalities remain valid, and it is the single deficit measure D0 that is summed. Conversely, without support disjointness one cannot replace this proof by an unpaid sum of local deficits.

## 2. One global normalization and the fixed layer geometry

Write L=logn and choose

\[
d=7/25,\quad J=128,\quad\sigma=3/50,
\quad x=(1-\sigma)L=47L/50,\quad t=n/x,\quad\beta=1+1/t,
\]

\[
\epsilon_0=2^{-48},\quad h=\epsilon_0/(Ld),
\quad\zeta=n^{-3},\quad\kappa=\beta-\zeta,
\quad\delta_r=2^{-30}.
\]

Partition [18/25,24/25] into 128 equal intervals with endpoints

\[
q_j=18/25+(6/25)j/128,\qquad 0\le j\le128.
\]

Assign shared endpoints to one adjacent source interval, for example by half-open intervals with the last endpoint closed. Let Ij be that interval with lower/upper endpoints qlo,j,qhi,j. Targets in layer j have radius

\[
r\in[r_{{\rm lo},j},r_{{\rm hi},j}]
=[q_{{\rm lo},j}+d+\delta_r,
q_{{\rm hi},j}+d-\delta_r].
\]

These target radial intervals have positive gaps 2δr between adjacent layers. Source radii are ρ=tq with q∈Ij. For large n all targets have 1<r<2 and all sources have ρ>4. Thus all Aj/Pj pairs are mutually disjoint in point space, with one common final support radius normalized to 1. There is no separate rescaling of individual layers.

The ambient corner measure is

\[
d\lambda(z)=2^{-n}m(z)^{-n}\mathbf1_{\mathcal D}(z)dz,
\quad m(z)=\max(1,\|z\|_\infty),
\quad\mathcal D=\{\|z\|_\infty<t:\exists i,z_i<1-\beta m(z)\}.
\]

The full-source measure before restriction to D has each maximal face at radius ρ weighted by dρ/(2ρ), with independent Uniform[−ρ,ρ] free coordinates.

## 3. Local inverse profiles, active R, and the typical targets

Use matching |q−(r−d)|<h, and put

\[
R(q)=q+d-h,\quad a_0(q)=(\kappa-1)\rho=(1-\zeta t)q,
\quad g=R-a_0=d-h+\zeta\rho.
\]

For a target radius r the upper matching endpoint qb=r−d+h has R=r, so it is the exact-mean baseline for J. The two count endpoints are q±=r−d±h.

When h<δr/4, every matched q stays strictly inside its own source interval and

\[
R\ge q_{{\rm lo},j}+d+\delta_r-2h
\ge R_j:=q_{{\rm lo},j}+d\ge1. \tag{2}
\]

This active-R guard is essential in the first layer: the full unmatched source edge can have R<1. No local lemma requiring R≥1 is applied there. Such points have zero target column because no r in the trimmed target band can match them.

On the full source interval define

\[
\ell_R(z_i)=\min(1,(R+\kappa\rho-|z_i|)/(2R)),
\quad G(z)=\prod_i\ell_R(z_i)>0,
\quad b_j(z)=\min(1,e^{-\chi_jx}/G(z)).
\]

For large n R>1/2 even on unmatched edges, so this bounded extension is well defined. It is made before cap restriction. Set b_j=0 outside its full source band.

Use errors

\[
e=\epsilon_J=\epsilon_N=\epsilon_S=\epsilon_T=\epsilon_D=2^{-20},
\quad\epsilon_{\rm Ret}=1/4096,\quad\epsilon_A=1/65536.
\]

For clarity the one-coordinate inverse normalizer and the count variables are

\[
J(v;q)=\frac1{2\rho}\int_{[-\rho,\rho]\cap(v-\kappa\rho,v+\kappa\rho)}
\ell_R(z)^{-1}dz,
\quad \ell_*(q)=\frac{R+a_0}{2R}=1-\frac g{2R},
\]

\[
N_i(q)=1+\mathbf1_{\{|w_i|<a_0(q)\}},\quad T_N(q)=\sum_iN_i(q),
\quad E_\pm=n+(n-1)a_0(q_\pm)/r.
\]

Let Aj be the target annulus intersected with D and the three explicit events

\[
|\sum_i\log J(w_i;q_b)|\le2^{-20},\qquad
|T_N(q_\pm)-E_\pm|\le en/4.
\]

On each maximal target face the pin has Ni=1, and the n−1 other indicators are independent Bernoulli with parameter a0/r; hence each count variance is at most n/4. These are exactly the count laws used in the accepted R056 source, not a law imposed on actual μ. The local one-coordinate lemma applies on active matches: d≤1/3, q≥18/25>1/2, R≥1, g≤min(d,R/2), ρ≥4r, ζt≤h/2 eventually. It gives conditional log mean/variance O(nd³/t²)→0 and the a.e./AC q-derivative bound 64[d²+1tail]/ρ.

Here Kw≤n≤4nd is automatic because 4d>1. Therefore

\[
|\sum\log J(w_i;q)|\le2^{-20}+1280\epsilon_0,
\quad P_0=1-2^{-20}-1280\epsilon_0>0,
\]

\[
P_0\le\prod_iJ(w_i;q)\le1/P_0.
\]

The two count endpoint events have failure O(1/(e²n)), and the count/pin normalization calibration follows once 3h+1/n≤e/2. Both this and the log baseline hold with probability 1−o(1), uniformly over the finite layer system. The target cap complement is exponentially small on each target face. Consequently, if Lj is any rigorous lower bound for log(rhi,j/rlo,j),

\[
\lambda(A_j)\ge m_jn,
\quad m_j=(1-\epsilon_A)L_j. \tag{3}
\]

There is no union probability loss across layers: (3) is a separate mass statement on disjoint radial sets, and the eventual threshold can be the maximum of finitely many thresholds.

Pointwise J→1 uniformly over all target-coordinate values and matched q. Hence eventually J∈[1−e,1+e]. The exact inverse-captured source-shell coefficient is

\[
U(w,\rho)=\frac{\prod_iJ(w_i;q)}{2\ell_*(q)}
\sum_i\frac{N_i(q)}{J(w_i;q)}.
\]

The preceding endpoint events and monotonicity give
\(1-e\le T_N(q)/(2n\ell_*)\le1+e\): the calibration errors are bounded by 3h+1/n, and the endpoint concentration adds at most e/4. This follows directly from R=q+d−h, R(q−)=r−2h, R(q+)=r, and \(\ell_*'=(1-\zeta t)(d-h)/(2R^2)>0\). Combining this count ratio with the product and individual-J bounds gives the same constants in every layer:

\[
U_{\rm lo}=P_0(1-e)/(1+e),\quad
U_{\rm hi}=(1+e)/[P_0(1-e)],
\quad U_{\rm lo}n\le U(w,\rho)\le U_{\rm hi}n. \tag{4}
\]

All these assertions follow from local coordinate estimates and full-face target laws; no distribution assumption is made about the actual μ.

## 4. Exact layer exponents and per-face inverse retention

Let phibar,j be a rigorous upper bound for

\[
\phi_j=-\log(1-d/(2R_j)),\qquad R_j=q_{{\rm lo},j}+d,
\]

and let mubar,j be a rigorous upper bound for

\[
\frac{R_j}{q_{{\rm lo},j}(1-e)}\,\overline\phi_j^{\,2}.
\]

Take a fixed rational χj>mubar,j, and define θj=exp(−χjx), kj=nhθj. On every active match, (2) and g≤d imply Qpin≤φj≤phibar,j. Also

\[
R/q=1+(d-h)/q\le R_j/q_{{\rm lo},j}.
\]

Thus on each captured source face, using the full inverse-capture law and independent free coordinates,

\[
E(V\mid\mathrm{face})\le\overline\phi_j+x\overline\mu_j,
\qquad\operatorname{Var}(V\mid\mathrm{face})\le xV_j,
\]

\[
V_j=\frac{2R_j}{3q_{{\rm lo},j}(1-e)}\overline\phi_j^{\,3},
\qquad V=-\log G.
\]

These follow from the exact free-coordinate integrals Efull[φ/ell]=RQpin²/ρ and Efull[φ²/ell]=2RQpin³/(3ρ), divided by J≥1−e. The pin is deterministic. No variance of a face mixture is equated to a within-face variance.

Put Δj=χj−mubar,j>0. For example the explicit eventual conditions

\[
x\ge2\overline\phi_j/\Delta_j,
\quad x\ge4V_j/(\epsilon_{\rm Ret}\Delta_j^2)
\]

give, on every captured face, P(V>χjx)≤εRet. Mixing those probability inequalities preserves the same bound. Soft rows before cap deletion consequently have shell coefficients between (1−εRet)Ulo nθj and Uhi nθj. Every target matching interval has exact ordinary width 2h within its source band, so its log-width is between 2h/qhi,j and 2h/qlo,j.

## 5. Layer-specific cap loss, τ, and local degrees

Require the two strict rate conditions

\[
\eta_j:=\frac{1-q_{{\rm hi},j}}{2q_{{\rm hi},j}}>\chi_j,
\qquad \sigma/2-(1-\sigma)\chi_j>0. \tag{5}
\]

The original corner-captured source coefficient is at most (n/2)exp[−(1−1/n)ηjx]. Since b_j≤1, removing that cap part costs at most its unweighted measure. Dividing by kj leaves a fixed multiple of h^-1 exp[−((1−1/n)ηj−χj)x], which tends to zero by (5). Hence the loss is at most εS kj eventually in every layer.

Choose the one global threshold τ=√[A0(B0−γn)], with A0=β^n and a fixed desired γ>0 below the total certified deficit. The original scalar estimates give B0∼nL and B0≤nL, so this is positive eventually, and

\[
\tau/k_j\le\frac{\sqrt L}{h}
\exp\{-[\sigma/2-(1-\sigma)\chi_j]L\}
\longrightarrow0.
\]

Thus τ≤εT kj eventually for every j. Because J is fixed finite, the maximum of all geometric, retention, cap and τ thresholds is finite. This is not a claim of a numerical threshold or validity at R056's old L0.

The refined target-column inequality Ni/(2ell_i)≤1 gives the weighted column≤2kj whenever it can be nonzero; then R≥Rj≥1 by (2). Inactive source edges have column zero. Row upper degree is at most 2Uhi kj/qlo,j. Consequently it is valid to take

\[
M_j=2U_{\rm hi}/q_{{\rm lo},j}>2.
\]

Here the actual final Borel kernel, to remove any support ambiguity, is

\[
P_j=\mathcal D\cap\{z:\|z\|_\infty/t\in I_j\},\qquad
K_j(w,z)=\mathbf1_{A_j}(w)\mathbf1_{P_j}(z)
\mathbf1_{\{|\rho/t-(r-d)|<h\}}
\mathbf1_{\{\|z-w\|_\infty<\kappa\rho\}}b_j(z),
\]

where r=||w||∞ and ρ=||z||∞. Its row and column degrees are respectively ∫Kj(w,z)dλ(z) and ∫Kj(w,z)dλ(w). The inverse-capture calculation was on the full source band, with b defined there before deletion of its cap complement. The ambient row lower degree minus the same global τ is at least bj kj, where

\[
b_j=(1-\epsilon_{\rm Ret})2U_{\rm lo}/q_{{\rm hi},j}-\epsilon_S-\epsilon_T.
\]

If bj>0, (1), (3) give local deficit at least mj bj n/Mj. Any exact number

\[
0<\gamma_j\le(1-\epsilon_D)m_jb_j/M_j \tag{6}
\]

is therefore an allowed conservative layer budget.

## 6. Actual finite histories, cross-band increments, and finite-a scaling

Use the original dilation mathfrak a≥8n^4, α=mathfrak a/(mathfrak a+1), and ccap=(mathfrak aβ−1)/(mathfrak a+1). At a fixed load support point normalize the last minimum positive-contribution radius to 1, and reflect its center into the positive orthant. The accepted actual uncaptured-increment measure is dominated by the finite corner envelope. Restrict both measures to the limiting cap D; write these as μa and λa. The full-envelope missing mass is at least λa−μa on D. Since the finite corner threshold exceeds 1−βm by [((β+1)m−1)/(mathfrak a+1)]>0, the limiting cap is contained in that envelope and λa=α^-nλ exactly on all Aj/Pj. All targets have ||w||<2 and every actual target item meeting them has radius<3. For a source increment point z=y/(mathfrak a+1), the common support 0∈mathfrak a Qi and y∈Qi give ||z||∞<ri. Thus its actual item radius is at least ρ>4. These target/source radius facts and the uncaptured-center bound are the same actual-increment statements used in R056 §8; none depends on which radial layer contains the point. The guard κ+2/(mathfrak a+1)<αβ therefore makes every included source precede and capture the same actual target center. The bounded weight bj(z)≤1 preserves its one original eligibility budget τ.

If an original increment crosses two source bands, each local μa(Pj) integrates only that part of the increment. Summing disjoint parts never exceeds its original wi. There is no assignment of full wi to every layer and no requirement that a small cube stay inside one bin. The same pointwise measure argument handles a target increment meeting multiple target bands.

Ambient energies scale by ca², degrees and target masses by ca=α^-n≥1. More explicitly the finite-a local deficit is at least
\(\lambda(A_j)[c_a r_j-\tau/k_j]/M_j\ge\lambda(A_j)[r_j-\tau/k_j]/M_j\).
Thus the same lower bound is valid at finite a; no unproved convergence of μa or selectors is used. Apply the direct-sum lemma to the disjoint sets Bj=Aj∪Pj; source shared endpoints are null for λa and μa. The total missing cap mass is at least

\[
n\sum_j\gamma_j
\]

uniformly for every mathfrak a≥8n^4 once n exceeds the fixed common threshold. These are inequalities on the same actual μa and the same history; a family of independent worst positions is never presumed simultaneously attainable.

If γ<Σγj, the original load, cover and common-f witness argument gives

\[
C_n\le\left(1+\frac{(\beta+1/\mathfrak a)^n}{\tau}\right)
\left[\tau+\mathfrak a^{-n}+\alpha^{-n}J_{\mathfrak a}(t)-\gamma n\right].
\]

Let mathfrak a→∞ only in these uniform scalar bounds. This yields

\[
C_n\le F(A_0,B_0-\gamma n). \tag{7}
\]

It retains general L1 inputs, every level and every positive common radius. No scale-truncated operator or changed witness family is substituted.

## 7. One fixed exact 128-layer construction and its complete cost

The record `upper/round062/compute/candidate_certificate.json`, SHA256 `b889e35ef909bd979937752c67410bf3913eb76cdf7dc6c720556e093d4c0711`, supplies the following exact construction. Each log is enclosed with 48 positive atanh-series terms and the remainder bound 2y^97/[97(1−y²)], then rounded outward to denominator 2^200. Define the rational phibar,j and mubar,j by those upper bounds, and set

\[
\chi_j=10^{-6}\left\lceil10^6\overline\mu_j\right\rceil+10^{-4}.
\]

These are exact record fractions, not displayed decimal approximations. Take Lj as the corresponding downward-rounded lower logarithm of rhi,j/rlo,j, and gammaj as the downward-rounded value of the right side of (6). The independent author verifier reconstructs every fixed row; it does not import the producer helper or replay its exploratory grid.

The exact verified assertions are: all bj>0; χj>mubar,j; both rate inequalities (5) hold in all 128 rows; and

\[
\sum_j\gamma_j>
\gamma:=429/2000.
\]

Therefore the geometric theorem, if independently accepted, makes (7) an original-operator bound at this γ. The complete scalar penalty is

\[
-\log(1-\sigma)\le\frac{\sigma}{1-\sigma}=3/47,
\]

\[
\gamma-3/47=14163/94000,
\qquad\gamma-3/47-3/20=63/94000>0.
\]

Also A0≤ne^-σL, B0≤nL and their extra positive root terms, divided by n, are at most e^-σL+2√L e^-σL/2+e^-L→0. Thus, for all sufficiently large n,

\[
\boxed{\quad
C_n\le F(A_0,B_0-(429/2000)n)
<n\log n-n\log\log n-\frac{3n}{20}.\quad} \tag{8}
\]

This is a proposed full operator conclusion, not a claim from the conditional record alone. It has no asserted finite dimension threshold and no global parameter optimum. At dimensions where only earlier explicit theorems apply, retain those bounds; take the minimum of all applicable results.

## 8. Evidence, scope and unresolved quantitative task

`verify.py` performs a new exact Fraction reconstruction of this single fixed 128-layer certificate, including all log remainder directions, rounded moment bounds, record fractions, row budgets, strict rate gaps and the clean γ gate. This is not a production search or replay of any old parameter grid. Its input identity and results are in `verification.json`.

The analytic result uses a fixed finite J. Making J or a record depend on n without uniform analysis is not covered. The proof also does not permit overlapping point supports to be summed freely. A future explicit threshold must bound the maximum of the fixed layer-retention and rate thresholds as well as the target probability guards; no such threshold is inferred from the R056 theorem or the finite scalar record.


The imported original-operator machinery is precisely the accepted R056 `upper/round056/refined_inverse/PROOF.md` §§2–3 (local d≤1/3 laws and face calibration), §4 (inverse shell and exact target-column factor), and §8 (actual increment, original finite-cover/witness bridge). Its frozen source identity is `6986ed9ec48b46f3136a9a56722d0092db4f2942d52e35c9f4de37eb279f399a`. The new arguments here are layerwise parameters, their active-R guard, per-layer retained energy, and summation against a single pointwise deficit measure.
