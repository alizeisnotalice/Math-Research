# R064 independent whole-operator review

Decision: **ACCEPT for the fixed-tuple eventual conclusions**, with the
accepted original finite-cover/corner-envelope/scalar interfaces retained.
For every sufficiently large integer n, the audited argument proves

`C_n <= F(A_n(t), B_n(t) - (61/200)n)` and
`C_n < n log n - n log log n - n/5`.

There is no certified explicit starting dimension. In particular no earlier
tuple's explicit threshold is transferred here. The endpoint coefficient
`gamma - sigma/(1-sigma)` is not claimed as a strict final improvement.
This review does not rely on finite numerical patterns establishing an
operator bound: the operator transfer is checked analytically below.

## Frozen identities and work performed

The complete author source actually read is
`/Users/zhengzhihao/Desktop/T/cube-r040-cover-selector/cube-dim-7-bounds/upper/round064/radial_calibration/PROOF.md`,
SHA256 `17ce38ab5af46c0358046bdcc2aa57a05844c0799e9a97280d959887f59615d6`.
The author's `certify.py` and `certificate.json` were read, but neither was
executed/imported as this audit's independent calculation. Their identities
are recorded in `source_lock.json`.

The new relative-coordinate lemma actually read in full is
`/Users/zhengzhihao/Desktop/T/cube-r030-roots-audit/cube-dim-7-bounds/audit/round064/relative_coordinate/LEMMA.md`,
SHA256 `c60d48d116863c9c9cac177e19a9fa4158f12fd8dcef0d5239a5ad0a5487d9f1`.
It is bound with its own active relative domain and laws; no old absolute-gap
lemma or old derivative sign is substituted. The original finite-history
cover/load construction and scalar definitions were also checked against
the private clone's `paper/main.tex` and
`revision_v3/asymptotic/proof.md`; these are accepted dependencies rather
than newly claimed improvements.

The new local Python run used exact Fraction and SymPy, with **41 checks
PASS**. Its target-mass log enclosure is a different recipe: 40-even/41-odd
alternating log1p sums, rather than the author's atanh series. It additionally
checks the new matching/count/measure and finite-dilation identities.
`verification.json` records actual interpreter, time, PID, script hash,
exact residuals and source identities. This run did not import an author
helper, rerun a grid or previous baseline, use shared MATLAB, or perform Git
operations. Only this audit directory was written.

## 1. New domains and genuine relative-coordinate binding

Here `R=(3/2)q-h` and `g=(1/2+epsilon)q-h`. The full source band is defined
before cap restriction and soft weighting, including inactive points.
The condition `epsilon qmax <= h/s` gives the exact identity and direction

`g - R/3 = epsilon q - 2h/3 <= 0`.

Together with `q>=2/3` and `h<=1/1000`, this gives `0<g<=R/3`.
It also forces `epsilon<=h`; all new derivative bounds are inside the
relative lemma's domain. An active source has `R>rlo-2h>1`; inactive
sources with `R<1` receive zero final column and are never inserted into
the local `R>=1` estimates. The full-band ell and soft weight are still
positive there because `R>0` and `a0>0`.

For each target, matching q ranges from `(r-h)/s` to `(r+h)/s` and stays
strictly inside the source band when `h<delta/2`. The right endpoint is
the mean-one baseline: `R(q_b)=r`. The reference radius at the left
endpoint is `r-2h`. The old additive baseline cannot be used.

The local lemma supplies the unconditional a.e. derivative/absolute
continuity bound with `R'=3/2` and `g'=1/2+epsilon`; no K_w event remains.
The product bridge over the whole q interval is at most
`512Lh=2^-55<2^-50`. At the baseline, the free target coordinates are
independent, mean log and variance are uniformly `O(L^2/n)+O(L/n)`.
Thus baseline probability tends to one uniformly over faces/radii.
The deterministic individual-J estimate tends to one separately and
uniformly; no product of n individual error factors or union bound is
being substituted for the product-log argument.

## 2. Target calibration, mass and complete inverse rows

The exact captured pinned-source count is
`N_i=1+1_{|w_i|<a0}` (boundary equality still contributes one).
The pinned target coordinate contributes one, and the other n-1 counts
are independent Bernoulli additions at either fixed endpoint. Therefore
each endpoint failure has probability at most `4/(eps_N^2 n)`.

The independently differentiated pin factor satisfies
`ell_pin'=-(1-epsilon)h/(2R^2)<0`. Since Nsum is nondecreasing,
`Nsum/(2n ell_pin)` is nondecreasing. At the lower endpoint its mean
deficit is exactly

`[2h a_-/(r(r-2h)) + a_-/(nr)]/[1+a_-/(r-2h)]`.

Because `R_-=r-2h>=1` and `a_-<r`, this is at most `2h+1/n`.
At the upper endpoint the mean ratio is at most one. Fluctuations from
the two count events contribute at most eps_N/4 to each relevant side.
This validates the stated whole-interval calibration when
`2h+1/n<=3eps_N/4`, with the new endpoint directions.

The original cap failure on each full target face is at most
`(3/4)^(n-1)` eventually. Baseline, two count events and cap are the four
bad events. Their uniform total is eventually at most eps_A. Hence
`lambda(A)>=(1-eps_A)n log(rhi/rlo)>mL n` for this same good target set.

All signed source faces are included in the exact inverse coefficient
`U=(product J_i)/(2ell_pin) * sum_i N_i/J_i`. Product, count and
individual-J bounds therefore give `Ulo n<=U<=Uhi n` without dropping
pinned or omitted-coordinate terms. The source logarithmic width is
`log((r+h)/(r-h))`; s cancels only in this ratio. Its multiplication
by the bounded target factor `r/rhi` gives the two bounds in (19).

## 3. Soft retention and cap payment preserve the laws

On each fixed captured signed source face under the full-source
capture/inverse law, the pinned log ell term is deterministic and the
remaining coordinates are independent. The exact first and second
full-source integrals yield the stated mean upper rate below
`999/20000` and variance upper rate below `1/150`. These use `R/q<=s`,
not an old absolute bound on g. The pin cost is at most 1/5; for x>=8000
the remaining threshold gap is at least x/40000. Chebyshev is applied
on **each** face and its failure bound is then mixed. No mixture-variance
identity or equality of different face means is assumed.

The soft ratio relative to theta/G is `min(1,G/theta)`. It is exactly
one on `G>=theta`; hence at least `1-eps_Ret` of the theta-scaled
inverse row is retained. The target radial factor is constant in that
fixed-target source probability law and changes none of its moments.
This proves both full-source raw-row bounds in (24).

Cap deletion is paid afterwards, without a conditioned independence
claim. The wide source-band mass bound must use `Csrc=727/2000`.
For q<=qmax<1, the negative pinned source face lies in D and free cap
fractions give eta=91/1818. The strict rate gap is
`eta-chi=1/18180`. Both actual weights are <=1, so any row's cap loss
is bounded by unweighted cap-complement mass. After division by the
same `k=nh theta/rhi`, the ratio is
`(Csrc rhi/h) exp(-((1-1/n)eta-chi)x)`, which tends to zero.
The rhi factor is present. This proves the eps_S payment and (26).

## 4. Sharp column, same kernel and common degree budget

For every active source and matched target, `r>=R` and
`kappa rho-|z_i|>0`. Thus each target-coordinate capture fraction is
strictly above one half. A partial fraction has exactly one captured
target endpoint, and a full fraction has at most two. Consequently
`T_i/(2ell_i)<=1` for every coordinate, including the strict-capture
boundary upper-bound cases. Summing **all** target faces gives
`(1/2)sum_i T_i product_{j!=i}ell_j(r)<=n product_i ell_i(r)<=nG`.

The final inequality uses decreasing capture fractions for `r>=R`.
Now `(dr/r)*(r/rhi)=dr/rhi`, and the ordinary matched target interval
has width at most 2h (not 2h/s). Source soft weighting gives column
at most `2nh theta/rhi=2k`. Inactive columns vanish. Restrictions to
A and D only decrease the column. Since M>2, **this one final kernel**
has both row and column degrees <=Mk. Cover and load are not taken
from two incompatible kernels.

## 5. Original history, actual eligibility and finite dilation

The original accepted selection fixes increments and weights at their
selection times and has nonincreasing radii. At a positive load point,
normalizing the last minimum-radius contributing cube and reflecting its
center leaves actual increments dominated by the original corner
envelope. To make the corner domination explicit: after reflection u_i>=0,
an uncaptured point in `||z||<t` cannot satisfy the upper alternative
`z_i>u_i+c_cap m`, since c_cap>1. Its lower alternative is contained in
`z_i<alpha-c_cap m`; hence the common corner envelope dominates it.
On the chosen D, the exact positive threshold difference displayed
in Section 8 makes `lambda_a=c_a lambda`, `c_a=alpha^-n`.

At any actual target increment j in A, the original uncaptured-item
inequality and `||u||<alpha` give `r_j<3` eventually. An actual source
increment i in P has `r_i>=rho>4`. It therefore precedes j in the
**same** history. The checked displacement guard for every
`mathfrak a>=8n^4` gives genuine center capture in j's original
eligibility test. Included pieces contribute at most the original w_i;
both source soft weight and target radial weight are <=1. Thus the
actual weighted row is <=tau on actual target increments. The proof
does not substitute an arbitrary model measure for those increments.

For fixed-tuple eventual n, beta<2 and t=n/x give both
`mathfrak a>2t` and
`zeta-(beta+2)/(mathfrak a+1)>0` uniformly for every such dilation.
The energy identity with `nu=lambda_a-mu_a`, and disjoint A/P, bounds
the two missing-measure terms by one degree budget times total corner
deficit. Energy scales by c_a^2, degrees by c_a and target mass by c_a.
Since c_a>=1, and the row-minus-tau lower bound is positive, the
inequality in (30) has the claimed direction. It gives at least
Gamma_L n, strictly more than gamma n. This is a genuine improvement
of the original uncaptured load on every relevant original history.

## 6. Tau, full cover cost, L1 and all radii

The accepted integral estimate gives `B_n>=nL/2` eventually at
`x=(1-sigma)L`, so the chosen tau is positive. The ratio tau/k tends to
zero with the exact positive exponential gap `1/20000`; all remaining
factors are powers of L. This pays eps_T within the same row budget.

The original first-rejection cover factor and captured budget remain
in (31). The new deficit is inserted into its uncaptured term only.
All finite histories and every dilation >=8n^4 satisfy this numerical
bound. Sending dilation to infinity for fixed n therefore takes a
limit of **valid numerical upper bounds**, not a limit of selectors,
increments or a scale lattice. It yields (32) at the specified tau.

For arbitrary L1 input, use |f| and each strict-level-set point's
positive true witness radius. Its corresponding small open cube permits
a finite covering of every compact subset. Selected dilates retain the
true witness property. The history normalization works for every
positive terminal radius, and Lebesgue inner regularity gives the whole
strict level set. There is no compact-support, bounded-radius,
finite-model or common fixed scale restriction. Both weights staying
below one is essential to the genuine mu<=lambda/eligibility transfer.

The complete scalar F, not B alone, is expanded in Section 9. Exact
`c_raw=37249/181800` exceeds 1/5 by `889/181800`. The displayed positive
remainder tends to zero. It is eventually strictly smaller than this
fixed slack, proving the strict claimed final inequality. All finitely
many geometric, probability and budget guards are uniform and eventually
simultaneous; none of this argument supplies an explicit N.

## Independent exact margins and limitations

The actual new certificate gives:

- mass lower enclosure minus `1239/4000`: > `0.00000213367968334`;
- `gamma_budget-61/200`: > `0.00143821702573`;
- `Gamma_L` approximately `0.306737765625`, M approximately `2.00782778117`;
- mean-rate and variance-rate residuals strictly positive;
- cap gap `1/18180`, tau gap `1/20000`, final scalar slack `889/181800`.

The decimal values are displays of exact rational residuals in the receipt.
The author's mass floor `1239/4000` exceeds producer floor `3097/10000`
by `1/20000`; they are distinct conservative lower bounds, not identical
receipts. This audit independently certifies the author's floor and all
budget uses of it. No new production search or explicit threshold is
asserted. Final integration and research-ledger responsibility stays
with the coordinating author.

Reproduce only the new audit calculation with:

```sh
PYTHONDONTWRITEBYTECODE=1 /Users/zhengzhihao/.local/share/math-research-tools/venv/bin/python /Users/zhengzhihao/Desktop/T/cube-r030-soft/cube-dim-7-bounds/audit/round064/radial_calibration/verify_independent.py
```
