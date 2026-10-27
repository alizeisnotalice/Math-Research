# R069 independent complete-operator threshold audit

**FULL THRESHOLD ACCEPT**, locked to the R069 source
`d1e0420ba6ce509d2a1e9d5e0c6c01fd8aedd62c52adc3bc187b584fd646da10`.
For every integer n with the natural logarithm `L=log n>=2^30`, all
uniform guards of the accepted R066 mechanism hold simultaneously, and

`C_n <= Phi66 = F(A_n(t), B_n(t)-(317/1000)n) < T-(23/100)n`,

where `x=(917/1000)L`, `t=n/x` and `T=nL-nlogL`. This certifies a new
explicit validity threshold for the **same** complete Phi66 and operator.
It neither changes the R066 kernel/gamma nor retroactively changes any
R066 historical eventual/conditional receipt. The threshold is sufficient,
not claimed minimal. Main coordination retains final integration responsibility.

## 1. Actual read scope and frozen dependencies

The R069 threshold proof, its registration and verifier were read. The
author verifier was not executed or imported. The accepted R066 full proof
and full review, the relative-coordinate lemma, the new first-profile lemma
and its source binding, and the exact-price proof were actually read and
their interfaces checked. Final source identities are:

| Source | SHA256 |
|---|---|
| R069 threshold proof | `d1e0420ba6ce509d2a1e9d5e0c6c01fd8aedd62c52adc3bc187b584fd646da10` |
| R066 operator proof | `4b7b3210f063312527c08adab5b7f81dc4b2c0ecf879e466dc715eb56ea6ebd1` |
| R066 full operator review | `5a13e5a4364a2fbcc4b9d89fb2536ebeded2569033041017a7355237fd110ea8` |
| R064 relative-coordinate lemma | `c60d48d116863c9c9cac177e19a9fa4158f12fd8dcef0d5239a5ad0a5487d9f1` |
| R066 first-profile lemma | `67bdf01d9aa65c70bf058df08b527b904e954d19f3327ecbb43e5b5c316b5030` |
| R066 exact-price proof | `d058c87d90db9748ac05a6d42fd711f68546daa217eee5ef8106b9c6e1d689f0` |

The accepted fixed budgets `(1-epsA)log(rhi/rlo)>637/2000`, `M>2`,
`(1-epsD)Gamma>317/1000`, target mean `mu0<221/5000`, source mean
and variance constants `223/5000` and `1/150`, and exact price bound
`log(1000/917)<4333/50000` are reused with their accepted scopes.
The fixed/log programs and their prior baseline checks were not replayed.
Their actual old receipts were read for identity and constants, without
changing their status or hashes. Paths and hashes are in `source_lock.json`.

## 2. Different independent calculation and all-tail method

One new local Fraction run produced **46 checks PASS**; actual time, PID,
interpreter, code/source identities and exact residuals are in
`verification.json`. There is no floating inference. In particular, this
audit does not depend on the author's 70 green gates.

The independent compression deliberately uses the weaker global bound

`n=exp L >= L^8/8! >= N_*(L):=L^8/2^16`,

since `8!<2^16`, and also `n>=L^2/2`. It uses direct baseline
Chebyshev with denominator `(e-mean_bound)^2`, **Taylor12** for both
source-cap and tau payments, and **Taylor5** for the full scalar error.
The author's recipes instead use `L^8/8!`, Taylor10 and Taylor4.
Neither recipe is a sampled extrapolation: `exp y>=y^j/j!` holds for
all positive y, and the resulting powers are explicitly decreasing.

Here `L0=2^30` and `sqrt(L0)=2^15` exactly. For every `L>=L0`,
`N_*`, `N_*/L` increase, while h decreases as `L^-1`. The upper
epsilon is `1/(xN_*^2)`, so its ratio to h decreases as `L^-16`.
All claims below are therefore certified on the infinite tail, not only
at a single tested dimension.

## 3. Geometry and correct relative local domain

At L0 the independent weaker bound already gives `N_*/L>16`,
`qmin N_*/L>4rhi>4`, `N_*>=400`, `N_*>=65`,
`h<min(delta/2,1/1000)` and `epsilon qmax<2h/3`.
They persist as just described. Since `x<=L`, `t>=N_*/L`, and
since `x>1`, `t<n`. Thus beta<2, kappa>1, both annuli lie below t,
and target/source supports are disjoint.

All full-band references have `R=sq-h>0`, `a0=(1-epsilon)q>0`,
and `g>0`. Their weights are defined before cap restriction.
Only active sources invoke the local R>=1 estimates: for them
`R>rlo-2h>1`. Inactive references may be below one and their final
columns are zero. Both matching trims hold on the entire target band,
and `q_b=(r+h)/s` gives `R=r` exactly.

The retained positive epsilon term gives
`g-R/3=epsilon q-2h/3<=0`. Hence the actual domain is
`q>=2/3`, `1<R<=r<3/2`, `rho>=4r`, `0<g<=R/3<1/2`,
with `epsilon q<=h/s` and `epsilon<=h`. These are the new relative
lemma's genuine hypotheses at qmax=.917; no absolute g<=1/3 or old
qmax=.909 bound is borrowed. No K_w event is needed because the
accepted log-J derivative is unconditional in the target coordinate.

## 4. Five good-target probabilities, including new Jphi

All target probabilities are taken separately on each fixed full signed
face and each r, with n-1 independent free coordinates and one fixed pin.
They are not laws assumed for an unknown actual history measure.

For the baseline log sum, the relative moments and pin bound give

`M_*(L)=[(3/4)L^2/qmin^2+(1/2)L/qmin]/N_*`,
`V_*(L)=(3/2)L^2/(qmin^2 N_*)`.

The run verifies `M_*<e/2` and `V_* /(e-M_*)^2 < epsA/8` at L0.
M_* decreases as a positive sum of `L^-6,L^-7`; V_* decreases as
`L^-6`, while the positive Chebyshev denominator increases. The
baseline failure is therefore strictly below epsA/8 throughout.

Each endpoint count has failure at most `4/(e^2 N_*)<epsA/8`.
The new `ell_pin'<0` and increasing Nsum give an increasing ratio
`Nsum/(2n ell_pin)`. The lower-end shortfall is at most `2h+1/n`
and upper-end mean is at most one. The independent drift check
`2h+1/N_*<3e/4` makes the correct whole-interval count calibration
valid. The n-1 pin correction and endpoint directions are retained.

For the original target cap, `1/rhi-L/N_*>1/2` persists since
`L/N_*` decreases. Each positive pinned-face complement consequently
has every free fraction below 3/4; negative pinned-face complements
are empty. Since n>=65 and `(3/4)^64<epsA/8`, this is a **per-face**
bound, with no unearned factor one half from averaging signs.

The new profile uses its baseline **target** Fubini mean, bounded by
`mu0 x`, with `mu0<221/5000`. That direction requires the accepted
lower bound on log(6/5); it is not the inverse-source mean or an odd
upper-log substitution. The separate pin is at most `1/rho_b`.
The check `N_* qmin>=200000` gives a pin cost at most `x/200000`.
The event threshold `(221/5000+1/100000)x` is therefore above the
total mean by at least `x/200000`. Bounded Jphi and within-face
independence give variance at most `x^2/(n qmin^2)`, hence failure
at most `90000000000/N_*<epsA/8`. This decays as `L^-8`.

No mutual independence of the five events is required. Their union
cost is at most `5epsA/8<epsA`. Thus the same target set has the
accepted mass lower bound `lambda(A)>mL n` at every L in this domain.

## 5. Full-q transport, inverse source law and common degree budget

The accepted a.e./AC log-J bound integrates over the entire q interval
to `512Lh=2^-71<2^-60`. The baseline log tolerance therefore gives
the scalar P0 product interval. The **separate** uniform individual-J
bound `L/(3qmin N_*)<e` decreases as `L^-7`; it is not multiplied
n times as a replacement for the product-log event.

The Jphi AC bound gives transport at most `16xh`; the checked
`h<1/1600000` pays the additional `x/100000`. Thus the same good
target has `sum Jphi(q)<=(2211/50000)x` at every matched q.
For each captured inverse-source face, the free phi mean is exactly
`Jphi/J_i`. This transfers the deterministic profile to **every**
source face before conditioning on cap or soft weight.

The accepted source mean `<1/5+(223/5000)x` and variance `<x/150`
combine with `x>=1000` and
`chi-223/5000-1/5000=1/5000`. At L0 the per-face Chebyshev bound is
exactly `1953125/11538530304<1/4096`, with positive residual about
`0.0000748708`; it decreases as `L^-1`. The uniform face probabilities
are mixed only after this calculation. No mixture variance is set equal
to within-face variance, and no cap-conditioned independence is used.

Soft retention and the complete pinned-source inverse coefficient
therefore retain the accepted Ulo/Uhi bounds. The radial relative-width
check `1/[1-(h/rlo)^2]<1+e` persists as h decreases. It uses the exact
`(r/rhi)log((r+h)/(r-h))` combination and preserves
`k=nh theta/rhi`. The column still sums **all** target faces with
`T_i/(2ell_i)<=1`, then uses the exact cancellation
`(dr/r)a_rad=dr/rhi` over ordinary target width at most 2h. Its bound
is 2k, and M>2 covers both degrees for this one final bounded kernel.

## 6. Finite-n cap loss, tau and full scalar positivity

The source band requires `Csrc=751/2000` and its actual k-relative
prefactor `Csrc*rhi/h`. The exact rate is `eta=83/1834`, with gap
`eta-chi=47/183400`. Since n>=N_*>=400, the explicitly checked
`eta-chi-eta/400>1/10000` yields a uniform safe exponent. Put
`cS=(917/1000)/10000`. Full-source cap deletion is bounded by

`S12(L)=12! Csrc rhi/[h_c cS^12 L^11]`.

At L0 this is below `3.87e-19<e`, and it decreases as `L^-11`.
Both actual weights are at most one, so this bounds the original
unweighted cap deletion after the full-source inverse-law argument.
It does not retain independence after deleting the cap complement.

The original scalar integral bound `I<=2t/(n-1)` and x in [L/2,L]
give `B>=T-4n/L`. For L>=16, `log L<L/4`; the checked
`L/4-4/L-1>0` and gamma<1 imply `B-gamma n>nL/2>0`.
This certifies the square root domain here, without an unquantified
eventual residue. Also x>4 gives B<=nL and A<=n exp(-sigma L).

For the same tau, the exact exponent is
`nu=sigma/2-chi(1-sigma)=47/200000>0`. The new full-tail bound is

`T12(L)=12! rhi/[h_c nu^12 L^(21/2)]`.

At L0 it is below `4.21e-19<e`, and decreases as `L^-21/2`.
Thus epsS and epsT payments hold at every L>=L0. The same target
mass, row excess b*k and degree M*k give the accepted strict
Gamma>gamma energy deficit.

For clarity, the author's stated Taylor10 bounds are also correctly
directed: `Csrc*rhi/h<2^110`, `cS L0>2^16`, and
`10!/2^50<e`; similarly the tau prefactor is `<2^126`,
`nu L0>2^17`, and `10!/2^44<e`. Their powers are respectively
`L^-9` and `L^-17/2`. The different stronger Taylor12 audit recipe
above makes the independent run distinct without changing any budget.

## 7. Actual finite histories, ca scaling, general L1 and every radius

All remaining guards in R066 Sections 8–9 now hold numerically in the
same explicit domain. Since t<n and n>1, every dilation `a>=8n^4`
satisfies a>2t and ccap>1. Beta<2 implies the exact displacement
guard is greater than `n^-3-4/(8n^4)>0`. No additional dimension
threshold is left in this part.

The accepted actual uncaptured increment measure remains dominated
by the finite corner envelope; on D this envelope is exactly
`lambda_a=ca lambda`, `ca=alpha^-n`. The explicit positive threshold
difference proves D's inclusion, and the total finite-cap deficit
dominates the missing mass restricted to D. Actual target items have
radius <3 since `alpha+rhi<3`, and source items have radius >4 from
their true increment points and common load support. They precede
those same targets in the one original radius-ordered history.
The positive displacement margin gives genuine center capture in
the original eligibility test. Included increment pieces cost at most
their original weights; both new factors are bounded by one.

Consequently actual weighted row is at most a_rad*tau<=tau. This
is not an ambient random-law claim about actual mu. Disjoint A/P
and the one-missing-measure identity give energy scale ca^2, degree
scale ca and target mass scale ca. Since ca>=1, the finite deficit
is at least the limiting bound >gamma n with the inequality direction
in R069 Section 7. The cover factor, captured tau and terminal a^-n
remain present.

Thus every original finite history and every dilation >=8n^4 satisfies
the strengthened numerical estimate. Every strict-level-set point of
arbitrary L1 input has its own true positive common-radius cube witness;
the small open cubes cover compact subsets, and selected dilates remain
true witnesses. Normalization works for every positive terminal radius.
Taking a limit only of the valid numerical bounds at fixed n, followed
by the accepted inner-regularity passage, yields Cn<=Phi66 for the
original Lebesgue centered cube operator with all positive radii.
There is no support assumption, radius cutoff/lattice, selector limit,
independent witness assignment or alternate operator hidden here.

## 8. Exact price, full root cost and old-Phi comparison scope

The accepted exact price bound `price<4333/50000` leaves the exact
margin `317/1000-4333/50000-23/100=17/50000`. A, the positive
cross term and 1/n are all retained in the remainder. Using the new
Taylor5 recipe gives the global upper bound

`R5(L)=2/L^2+5!/(sigma L)^5+2sqrt(L)5!/(sigma L/2)^5`.

At L0 this is below `1.74e-18<17/50000`; its three powers
`L^-2,L^-5,L^-9/2` decrease. This proves the strict full-Phi
23/100 inequality for the whole explicit domain, not only a scalar
coefficient comparison. No endpoint `gamma-price` claim is made.

The old Phi64 comparison properly retains its positive B buffer:
`B64>=T+(91/1000)n-4n/L`. Therefore on the old expression's own
valid domain, `Phi64>=T-[61/200-91/1000+4/L]n`. Already at
L=2^20 that bracket is below 9/40, which is below 23/100; it then
decreases. Hence the asserted full-cost ordering holds on the
intersection of the new explicit domain and the old mechanism's
true eventual validity domain. It supplies no explicit N64 and does
not promote 2^20 to an operator threshold. The historical old
source/receipts retain their original identities.

## Reproduction and limitations

```sh
PYTHONDONTWRITEBYTECODE=1 /Users/zhengzhihao/.local/share/math-research-tools/venv/bin/python /Users/zhengzhihao/Desktop/T/cube-r030-soft/cube-dim-7-bounds/audit/round069/profile_threshold/check_compressed.py
```

Only this audit directory was written. No author helper, old log sum,
old passed baseline, production grid, shared MATLAB, Git or n=7 work was
run. There is no unresolved guard or requested source correction in
this threshold audit. The proof and new receipt certify this fixed
threshold and same-Phi consequence, not minimality or any new kernel.
