# R064 independent local audit: radial target calibration, complete degrees, same-history payment

**Decision: LOCAL ACCEPT.** The target radial weight and bounded source-inverse weight give the stated complete row/column coefficients and the same-history deficit interface. No missing factor s or rhi, pinned-face omission, or reference-radius reversal was found. This accepts only the local interface proved below, conditional on the independently reviewed coordinate/product and soft-retention inputs. It does not accept the complete new C_n theorem, its tuple activation threshold, its probability/scalar chain, or a numerical n/5 conclusion.

Frozen source read in full: `/Users/zhengzhihao/Desktop/T/cube-r040-cover-selector/cube-dim-7-bounds/upper/round064/radial_calibration/PROOF.md`.

Verified SHA-256: `17ce38ab5af46c0358046bdcc2aa57a05844c0799e9a97280d959887f59615d6`.

The separate relative-coordinate proof is cited by the source with hash `c60d48d116863c9c9cac177e19a9fa4158f12fd8dcef0d5239a5ad0a5487d9f1`; it was **not** re-audited here. This audit does not replace that dependency or the separate full-source retention/scalar review. No author verification script, old grid, LP, MATLAB, or production integration was run.

## 1. Pointwise definitions and required active domain

Fix n and the tuple s=3/2. At a source radius rho=tq set

    R=sq-h,
    a0=(kappa-1)rho=(1-epsilon)q,
    g=R-a0,
    ell_R(z_i)=min(1,(R+kappa rho-|z_i|)/(2R)),
    G_R(z)=product_i ell_R(z_i),
    b(z)=min(1,theta/G_R(z)).

At target radius r, matching is `|r-sq|<h`. The target factor is `a(w)=r/rhi` on A0, zero off A0. Therefore 0<=a,b<=1 pointwise, independently of any actual measure.

For fixed target r the **entire** matching interval is

    q_-=(r-h)/s, q_+=(r+h)/s,

inside the full source q band under the frozen band's margin guard. It has ordinary q-width 2h/s. At q_+, R=r; at q_-, R=r-2h. For every actually matched source, R>rlo-2h>1. The full source band also contains some inactive sources with R<1, but these have zero final column and are not used in active R>=1 coordinate inputs. Defining b on them is harmless: R>0, a0>0, and `kappa rho-|z_i|>=a0>0`, so their G_R remains positive.

All row statements below use matched q, hence active R in `[r-2h,r]`, not an arbitrary full-band R. All fixed-source column statements keep R fixed while integrating r>=R. This distinction is essential and is respected by the frozen proof.

## 2. Complete fixed-source target-face coefficient

Fix one source point z, with rho=||z||_infinity, and put `a_i=kappa rho-|z_i|>0`. The active guard rho>=4r, kappa>1 ensures kappa rho>r. At target radius r, the exact free-target-coordinate capture fraction is

    ell_i(r)=min(1,(r+a_i)/(2r)).

Let T_i(r) count the captured pinned TARGET signs w_i=+-r. Because a_i>0:

* If ell_i(r)<1, exactly one target endpoint is lost, so T_i=1.
* If ell_i(r)=1, T_i<=2. At an equality/tangency with strict capture, T_i may be one; this only improves the upper bound.

It follows in every case that `T_i/(2ell_i)<=1`. A signed maximal target face has normalized shell mass `dr/(2r)`. Summing **all** 2n faces gives the coefficient per dr/r

    C_face(z,r)=(1/2) sum_i T_i(r) product_(l!=i) ell_l(r)
              =G_r(z) sum_i T_i(r)/(2ell_i(r))
              <=n G_r(z).                               (1)

For fixed z, every a_i is constant and positive, so ell_i(r) is nonincreasing in r. Since matched r>=R,

    G_r(z)<=G_R(z).                                      (2)

This is stronger than the older generic deleted-coordinate bound by 2n G. It is justified here by the positive a_i and the exact pinned-target count; it is not assumed from a sampled angular average.

Now multiply the full face measure by the target radial factor:

    (dr/r)*a(r)=dr/rhi.                                 (3)

The full matched target r interval, after any intersection with A0, has ordinary width at most 2h. Equations (1)--(3) imply

    C_full,radial(z)<=2nh G_R(z)/rhi.

The source weight is constant for this fixed z, and `b(z)G_R(z)<=theta`. Therefore

    C_lambda(z)<=2nh theta/rhi=2k,                       (4)
    k=nh theta/rhi.

The good-target and corner-cap restrictions only reduce this full column. Inactive sources have zero column. No unaccounted 1/ell_pin or factor two remains in (4); pinned TARGET faces have already been included in (1).

## 3. Complete source-face row coefficient and the pin-count calibration

Fix target w with r=||w||_infinity in A. Let

    J_i(q)=(2rho)^(-1) integral_[-rho,rho]
                1_(|z_i-w_i|<kappa rho)/ell_R(z_i) dz_i.

On a pinned SOURCE face z_i=+-rho, its inverse coordinate weight is `1/ell_pin`, with `ell_pin=(1+a0/R)/2`. The number of captured pinned source signs is

    N_i(q)=1+1_(|w_i|<a0(q)).

One sign is always captured on the active domain rho>=4r, kappa>1; the other is captured exactly under the displayed strict inequality. Thus the coefficient per d rho/rho=dq/q, before theta soft calibration and target radial weighting, is exactly

    U(w,rho)=(1/(2ell_pin)) sum_i N_i product_(l!=i) J_l
            =(product_i J_i)/(2ell_pin) sum_i N_i/J_i.    (5)

This includes every pinned SOURCE face and every omitted-coordinate term. The pinned source factor is not silently replaced by a free-coordinate J.

The frozen endpoint count calibration is valid. With epsilon constant at fixed n,

    ell_pin'(q)=-(1-epsilon)h/(2R^2)<0,

while Nsum(q)=sum_i N_i(q) is nondecreasing. Hence `Nsum/(2n ell_pin)` is nondecreasing on the whole matching interval, including its step boundaries.

For a full target face, the pinned target coordinate has N_i=1, because `a0(q)<=q_+=(r+h)/s<r`; the other n-1 coordinates give mean `1+a0/r`. Thus endpoint means are exactly

    E Nsum(q_+/-)=n+(n-1)a_+/-/r.

At q_+, R=r and the mean ratio is at most one. At q_-, R=r-2h>=1, and direct subtraction gives

    1-[1+(1-1/n)a_-/r]/[1+a_-/(r-2h)]
      =[2h a_-/(r(r-2h))+a_-/(nr)]/[1+a_-/(r-2h)]
      <=2h+1/n.                                        (6)

An endpoint fluctuation of at most `eps_N n/4` changes the ratio by at most eps_N/4, since ell_pin>=1/2. Therefore the stated guard `2h+1/n<=3eps_N/4` proves

    1-eps_N<=Nsum(q)/(2n ell_pin(q))<=1+eps_N             (7)

for every matched q. Only the endpoint count events and deterministic monotonicity are needed for this extension. No wrong additive-gap monotonicity direction is imported.

Now take as separate coordinate/product inputs

    1-eps_J<=J_i<=1+eps_J,
    P0<=product_i J_i<=1/P0.

Using eps_N=eps_J=e in (5)--(7) yields exactly

    Ulo n<=U<=Uhi n,
    Ulo=P0(1-e)/(1+e),
    Uhi=(1+e)/(P0(1-e)).                                (8)

These are conditional degree consequences of the independent coordinate/product inputs, not a re-audit of the relative J formula, its AC bridge, or its baseline concentration.

## 4. Exact radial row-width normalization

For the full matching interval,

    integral_q_-^q_+ dq/q=log((r+h)/(r-h)).               (9)

The s factor cancels exactly here. It would be incorrect to replace this by its ordinary q-width 2h/s without also dividing by q. Put u=h/r. The positive atanh series gives

    2u<=log((1+u)/(1-u))<=2u/(1-u^2).

Multiplication by a(w)=r/rhi therefore gives

    2h/rhi <= (r/rhi) log((r+h)/(r-h))
             <=[2h/rhi]/[1-(h/r)^2].                   (10)

The source's eventual guard that the last multiplier is <=1+e proves its upper width bound, and its weaker lower factor 1-e is valid. The reference radius R varies with q in the row integral, as it must; (8) is uniform over that whole matching interval, so no constant-R integration is being substituted.

Assume the separately audited inverse-face retention input: the soft coefficient before target radial weighting is at least `(1-eps_Ret)theta Ulo n` and at most `theta Uhi n` at each matched q. At fixed w the target weight a(w) is a constant, so multiplying by it does not change the conditional inverse-source-face retention probability. Integrating (9)--(10) proves

    2(1-eps_Ret)(1-e)Ulo k <= raw radial soft row
                              <=2(1+e)Uhi k=Mk.         (11)

The lower 1-e is conservative. In particular the row and the column use the **same** `k=nh theta/rhi`. There is no leftover s or rhi in their ratio. Since Uhi>1 and e>0, M=2(1+e)Uhi>2, so Mk controls both the upper row and column (4).

## 5. Source-cap subtraction normalization

Suppose the separate full-source cap estimate is

    Lcap<=Csrc n exp[-(1-1/n)eta x].

Since a,b<=1, deletion of this set costs at most Lcap in each target row. Relative to the actual k in (4) and (11), the ratio is

    Lcap/k <=(Csrc rhi/h)
                 exp[-((1-1/n)eta-chi)x].              (12)

Thus the extra factor rhi/h is necessary and is correctly retained in the frozen proof. There is no cancellation that would allow the cap error to be divided merely by nh theta. Given the separate uniform cap gate making (12)<=eps_S, the final lower row is

    R_lambda(w)>=[2(1-eps_Ret)(1-e)Ulo-eps_S]k.          (13)

The cap exponent/probability theorem itself is outside this audit; (12) audits its normalization and placement in the degree budget.

## 6. Same-history actual row and finite-a deficit

All kernel factors are nonnegative and at most one, in particular the new **pointwise** target factor a(w)=r/rhi on A. It depends on a target point, not on a whole original item label. Restricting an item's increment to included kernel points can only decrease its mass contribution. Therefore the real-history capture/precedence bridge, if valid for the included pairs, still implies

    integral K(w,z) d mu_a(z)<=a(w)tau<=tau

on mu_a-a.e. target point. It is permissible to use the coarser tau payment in the energy argument; no factor rhi or 1/a(w) must be paid. The source factor is `min(1,theta/G_R)`, not an unbounded inverse weight.

The frozen bridge's geometry is consistent with this interface. Targets have point radius r<rhi<3/2 and actual item radius <3; source points have rho>4 and actual item radius at least rho. Thus every included source item precedes the same actual target item. The stated displacement guard reduces to

    zeta>(beta+2)/(mathfrak a+1),

which follows at mathfrak a>=8n^4 once beta<2. The actual center displacement then gives strict capture at that target's own eligibility test. This audit does not certify the remaining tuple gates or reconstruct the original envelope proof; these are retained bridge dependencies.

For clarity the energy algebra is independent of item labels. Let E be the ambient energy of the final **radially weighted** kernel, D=Mk its maximum degree, and lambda(A) its unweighted ambient target mass. Let `nu_a=lambda_a-mu_a`, with `lambda_a=c lambda` on A and P, where `c=alpha^(-n)>=1`. The two-measure identity is

    E_lambda_a,lambda_a-E_mu_a,mu_a
      =E_nu_a,lambda_a+E_mu_a,nu_a.

Both terms are bounded by the scaled degrees cD, and A and P are disjoint. Consequently

    Delta_a >=[c^2 E-tau c lambda(A)]/[c Mk]
             =[cE-tau lambda(A)]/[Mk]
             >=[E-tau lambda(A)]/[Mk].                  (14)

Assuming the separate threshold gate `tau<=eps_T k`, equations (13)--(14) give

    Delta_a>=lambda(A)b/M,
    b=2(1-eps_Ret)(1-e)Ulo-eps_S-eps_T.                 (15)

If the separate good-target mass input is `lambda(A)>=mL n`, the resulting local deficit is at least `mL b n/M`. The mass in (15) is not `integral a d lambda`: the target weight was already included in the uniform lower row (13). Replacing it again by a weighted target mass would pay the radial factor twice. Conversely dropping rhi from k before deriving both degrees would be invalid. The frozen proof does neither.

Thus no missing rhi factor occurs in the claimed same-history interface. The activation of its positive rational deficit, the chosen gamma, full probability inputs and final C_n scalar bound remains the responsibility of the separate whole-chain review.

## 7. Actual new bounded exact checks

`verify.py` was actually run under Python 3.15.0 using only Fraction arithmetic. It performs seven new local checks:

1--4. Four source-coordinate configurations in dimensions 1--4. It directly intersects target-coordinate intervals, explicitly counts **both pinned target signs**, sums all target faces, and checks `C_face<=nG_r<=nG_R`, including saturation and an exact tangency. It also checks the soft-weighted `2nh theta/rhi` column bound.

5. Exact q-band/baseline identities and the rational atanh-series upper/lower width envelopes; no floating logarithm is used.

6. A deterministic endpoint-count example with n=4096 represented by two coordinate multiplicities, checking the two endpoint count events and their complete ratio calibration. This is a constant-size rational calculation, not a high-dimensional integration or a typicality claim. No J formula is evaluated.

7. Exact cancellation of k and rhi in the finite-a energy expression with the source's rational Ulo/Uhi/b/M definitions, plus the cap/k normalization factor in (12).

All seven passed. The small coordinate cases check the local formulas' hypotheses, not the source's full fixed n,t,L tuple or an original operator history. The count example does not claim that its coordinate-product event holds; the product input remains separate. No old verification, root/J diagnostic, parameter grid, LP, production batch, or MATLAB was run.

`verification.json` records the actual outputs and hashes. Reproduction: run `python3 verify.py` in this audit directory. Only this independent clone's `cube-dim-7-bounds/audit/round064/radial_degrees/` was written. No other source or main clone was edited, and no Git operation was performed. The parent retains whole-chain proof and integration responsibility.
