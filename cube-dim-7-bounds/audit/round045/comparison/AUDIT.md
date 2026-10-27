# R045 comparison and d-endpoint: independent audit

**Decision: ACCEPT for the conditional comparison and endpoint corollary actually stated.** No counterexample or blocking gap was found. This decision does not independently accept either new full-operator mechanism or assign either one an explicit operator threshold.

Frozen audited source: `/Users/zhengzhihao/Desktop/T/cube-r030-parameters/cube-dim-7-bounds/upper/round045/comparison/PROOF.md`.

SHA-256: `ba7498d1ccdc2c9158b065ade728c60f7b83cc2706ce517580d866903451d2f9`.

The source was read in full. The two mechanism sources were used only to identify their parameter/cost statements and d-dependent guards. Their full-operator proofs remain separate dependencies. Neither the author's `verify.py` nor its receipt was executed or reproduced. An independent compressed exact check was actually run, as described below. No shared MATLAB, production calculation, Git synchronization, author-source edit, or main-clone edit was performed.

## 1. Dependency and activation boundary

The geometry source has hash `add435fdf2f55993eba9064bd50c493e4c42769c1c48463329c4015be2e66ab9` and explicitly begins with **every fixed** `0<d<=2^-24`. Its proposed cost is

    F(A_n(n/[(1-32d^2)L]), B_n(n/[(1-32d^2)L])-nd/65536).

The soft source has hash `22db93b275416b90e9f15a0118a88cd412fee0c2d0b20c4cce1ed32769f259ac`, fixed `d=2^-64`, `B=4K=2^22`, and cost

    F(A_n(n/[(1-Bd^2)L]), B_n(n/[(1-Bd^2)L])-2^-32 nd).

Accordingly the comparison's I/S parameters are exact:

| Formula | sigma | Actual Delta | Display Delta/2 |
|---|---:|---:|---:|
| I at d=2^-30 | 2^-55 | 2^-46 | 2^-47 |
| S at d=2^-64 | 2^-106 | 2^-96 | 2^-97 |
| I endpoint at d=2^-24 | 2^-43 | 2^-40 | 2^-41 |

Scalar ordering does not require the new operator mechanisms. An inequality involving C_n additionally requires independent acceptance of the relevant full chain and entry into that chain's own valid domain. In particular neither `L>=2^240` nor old R042's explicit domain certifies the new mechanisms. For the endpoint one must accept the geometry theorem over its stated fixed-d range, not merely accept its displayed d=2^-30 specialization. The source preserves all these conditions.

## 2. Complete scalar envelopes

Write `T=nL-n log L`, `x=(1-sigma)L`, and `t=n/x`. For `x in [L/2,L]`, `L>=64`, the accepted cap identity and integral bound give

    B >= n log t - (n/2)*(2t/(n-1))
      >= n(L-log x)-2n/x >= T-4n/L.

The core term `1-h(t)^n` is nonnegative because t>1. The upper bound is `B<=1+n(L-log x)`, and x>e implies `B<=nL`. Also `A<=(1+ x/n)^n<=exp(x)=n exp(-sigma L)`.

The lower envelope for the **full** cost is therefore

    F(A,B-Delta n) >= B-Delta n >= T-n(4/L+Delta).

It is valid only when the root argument is nonnegative; the source supplies a stronger common condition: `B>8n` and every new Delta<1. For example `log L<=L/2` on this domain gives `B/n>=L/2-4/L>8`. Thus all new roots are positive. This does not discard a square-root term of unknown sign.

For the upper envelope expand the complete root cost:

    F(A,B-Delta n)=A+B-Delta n+2 sqrt(A(B-Delta n)).

Use `-log(1-sigma)<=2sigma`, `A<=n exp(-sigma L)`, and `B-Delta n<=B<=nL`. This yields exactly source (3), including the positive `2sigma n`, the constant 1, and both root/exponential terms. The comparison never replaces the actual Delta by its display coefficient Delta/2.

## 3. Explicit full-tail control and new-cost ordering

Taylor positivity gives `exp(u)>=u^2/2` for u>0. Applying this to L, sigma L and sigma L/2 gives

    1/n+exp(-sigma L)+2 sqrt(L) exp(-sigma L/2)
      <=2/L^2+2/(sigma^2 L^2)+16/(sigma^2 L^(3/2)) = E_sigma(L).

All three powers decrease on the entire positive L axis. The source's exact dyadic exponents at `L0=2^240` are correct. An independent, coarser certificate is enough:

| Formula | A strict common upper bound for E_sigma(L0) | Available margin |
|---|---:|---:|
| I | `3*2^-246` | `2^-47-2^-54` |
| S | `3*2^-144` | `2^-97-2^-105` |
| I endpoint | `3*2^-270` | `2^-41-2^-42` |

Each listed upper bound is strictly below the corresponding margin, as checked by Fraction arithmetic. Decreasing E_sigma extends the conclusion to **all** L>=L0, rather than just the base point. Hence the full-cost inequalities `Phi_I<T-2^-47 n` and `Phi_S<T-2^-97 n` are proved.

To order I and S themselves, use S's **actual** Delta=2^-96 in its lower envelope:

    Phi_S >= T-n(4/L+2^-96).

For L>=L0, `4/L+2^-96 <=2^-238+2^-96<2^-95<2^-47`. Combining this with I's full upper envelope proves `Phi_I<Phi_S`. This is a complete-cost comparison across two different x values; it does not compare only two displayed upper bounds.

## 4. Old complete costs and their individual domains

R042 uses `x42=L-L^(1/3)`, its own `t42=n/x42`, and its **actual** budget `2^-31 n L^(-1/3)`. Its lower envelope is

    Phi42 >= T-n[4/L+2^-31 L^(-1/3)].

At L0 the bracket is `2^-238+2^-111<2^-110<2^-97`, and both terms decrease. The x42 range required by the cap bound holds here. The source also correctly retains R042's already accepted explicit operator condition `L>=2^108`, which is covered by L0. Thus `Phi_S<Phi42`, and I is smaller still.

For the older listed fixed-t instances, `x_old=L-3 log L` lies in `[L/2,L]` throughout this comparison domain. Their actual budgets satisfy

    D_old <= 2n sqrt(log L/L)

on their respective valid domains: this uses 1/128 for R033, the theoretical `q<q_inf<1/8` for R034, `2^-16` for R035/36, and the actual integrated-energy mechanism coefficient `<sqrt(2)<2` for R040. It does not use a display coefficient in place of that last mechanism coefficient. R034's theoretical q bound already handles all its legal instances; the precise maximum of twelve recorded examples is not needed to infer the common lower envelope.

Consequently all these complete costs obey

    Phi_old >= T-n[4/L+2 sqrt(log L/L)].

The square-root term decreases when log L>1. Since `log L0=240 log 2<240<256`, its base value is strictly less than `2^-115`; the full bracket is less than `2^-114<2^-97`. This proves the claimed strict ordering throughout the tail.

The source correctly keeps each old instance's independent validity conditions. The elementary `log 2>2/3` yields `log L0>160>128`, covering R036's explicit condition. The recorded R034 Y values are 14 through 256, so their Y^2 conditions are covered; an arbitrary theoretical R034 instance still needs its own Y^2 condition. R035's original eventual threshold and the fixed-parameter R040 thresholds are not silently replaced by L0.

These lower envelopes are lower bounds for **formula costs**, not lower bounds for C_n. Taking a minimum over applicable old bounds and preserving their domains remains legitimate.

## 5. Unit-core all-t comparison

This is deliberately a separate eventual argument. Under the already accepted old chain, `V_unit>=U_old-2` and

    U_old=T+(3 log L+1)n/L+O(n(log L)^2/L^2).

The remainder divided by the positive correction is `O(log L/L)->0`, while that positive correction tends to infinity. Thus eventually `U_old-2>T`, and therefore `V_unit>T`. Both new scalar formulas are below T on the explicit comparison tail. Their new operator own-N conditions, and the old asymptotic threshold, remain additional separate requirements.

No bound for the old O constant is supplied, so L0 cannot be advertised as an explicit all-t unit-core threshold. The source does not make that advertisement.

## 6. Fixed-d guard audit, including the closed upper endpoint

The following checks address inheritance of the accepted parameter **range**, not acceptance of the inverse-row/operator chain itself. Set `dmax=2^-24`, `epsilon0=2^-20`, `sigma=32d^2`, and fix any d in `(0,dmax]` before sending n to infinity.

| Guard group | Reason it remains valid at every fixed allowed d |
|---|---|
| Scalar fraction | `0<sigma<=2^-43<1/2`, so x is between L/2 and L and t>1. |
| Matching widths | `h/d=epsilon0/(Ld^2)->0`. Also `zeta t/h=d/[epsilon0(1-sigma)n^2]->0`. Thus eventually h<d/16 and zeta t<h/2. |
| R, g and q | `q>=1-3d/4>1/2`; `R>=1+d/4-h>1`; `g=d-h+zeta tq` is positive and <=d. With R>1 and d<=dmax, `g<=R/2`. Also g'=zeta t<=d eventually. |
| Source-target separation | The target r is bounded by `1+3d/4<3`; the source rho=tq tends to infinity, so rho>=4r (and rho>=4R) and source radius>4 eventually. `kappa>1` follows from zeta t<1. |
| Coordinate-law conditions | From g<=R/2 and rho>=4R, `g/rho<=1/8`. Matching gives `T<=g+2h<2d`. The frozen coordinate bounds therefore retain their own stated local hypotheses. |
| Typical targets | The free-coordinate log mean and variance are `O(nd^3/t^2)=O_d(L^2/n)->0`; the pinned term tends to zero. Also nd->infinity, so the Kw threshold 4nd dominates `1+2nd`. These are fixed-d statements, not uniform assertions as d tends to zero. |
| Target cap | A free target coordinate misses D only with fraction `1/2+u/(2r)+1/(2t)<3/4` eventually, since u<=3d/4. Thus the excluded corner part remains exponentially small in n at the endpoint. |
| Matching bridge | `1280 L d h=1280 epsilon0=5/4096<1/8`, independent of d. The bound uses n/rho<=2L and d<=1, both retained. |
| Inverse-row Markov | `d+2Ld^2<=3Ld^2` requires only Ld>=1, eventually true for every fixed positive d. |
| Source cap deletion | `12d<=12dmax<1/64`; hence `d/32-12d^2>d/64>0`. The cap/k bound is `(Ld^2/epsilon0) exp[-L(d/32-12d^2)]`, which tends to zero. |
| Actual eligibility and tau | b<=1 remains true by definition. The computed ratio `tau/k <= (d/epsilon0)L^(3/2) exp(-4d^2 L)` tends to zero for each fixed d>0, so it eventually meets the required 3/128 threshold. |
| Finite a | With beta<2 and a>=8n^4, `zeta(a+1)>beta+2` eventually (in fact the displayed lower bound is already 8n). This is exactly `kappa+2/(a+1)<alpha beta`. Also a>2t eventually. These guards introduce no special d<=2^-30 requirement. |
| Deficit margin | The mechanism's lower deficit is `3nd/65536`; the asserted subtraction is only `nd/65536`. The coefficient gap is independent of d. Its finite-a scaling argument is inherited from the conditionally accepted mechanism, with no endpoint-specific sign reversal. |

No displayed condition needs d<=2^-30. The endpoint is included in the frozen hypothesis d<=2^-24; it is not a limiting passage from open d values. All these guards may require a **different** sufficiently large n when d is changed. This audit does not identify that n with the one for d=2^-30 or with L0.

## 7. Endpoint complete cost and correct parameter substitution

The endpoint recomputes the **entire** tuple `(d,sigma,x,t,h,theta,k,delta,tau)`. It does not keep the old d=2^-30 scalar parameter while merely increasing its deficit. At dmax,

    Delta=d/65536=2^-40,
    sigma=32d^2=2^-43,
    2sigma=64d^2=Delta/4.

Thus the scalar positive buffer uses only Delta/4. To obtain a Delta/2 display saving, there is still Delta/4 left for the constant and exponential/root terms. These remaining normalized terms are o(1) at each fixed d. The full root cost therefore eventually satisfies

    Phi_I,max <= T-2^-41 n.

On L>=L0, the independent dyadic envelope of Section 3 proves the stronger strict scalar inequality as well. Root positivity remains available since Delta<1 and B>8n. Thus the endpoint's budget B-delta and the optimized tau are legitimate.

For the old displayed I instance, the lower envelope at **its own** x gives `Phi_I>=T-n(4/L+2^-46)`. On L>=L0, `4/L+2^-46<2^-45<2^-41`, proving the endpoint full cost is strictly smaller. Both its improved C_n implication and this operator-level superiority activate only inside the relevant independently accepted own-N domains.

The source does not claim that this d is globally optimal, or allow d to vary with n while retaining a threshold. Those restrictions are necessary and are preserved.

## 8. Actual independent computation and source hashes

`microcheck.py` was actually run under Python 3.15.0, using only Fraction arithmetic. It checks the three **coarser dominant-term** envelopes of Section 3, the common strict comparison gaps, the endpoint equality `2sigma=Delta/4`, the cap exponent and bridge constants, and simple algebraic guards at a different probe `L=2^256` for d=2^-32, 2^-30 and 2^-24. All checks passed. The probe is explicitly not a certified full-operator threshold. No exponential evaluation, cap integration, old production program, parameter search, or author verification/receipt was run.

Reproduce by running `python3 microcheck.py` in this audit directory. `microcheck.json` records the exact rational values and the hashes below:

| Read dependency | SHA-256 |
|---|---|
| New comparison | `ba7498d1ccdc2c9158b065ade728c60f7b83cc2706ce517580d866903451d2f9` |
| Geometry mechanism (conditional) | `add435fdf2f55993eba9064bd50c493e4c42769c1c48463329c4015be2e66ab9` |
| Soft mechanism (conditional) | `22db93b275416b90e9f15a0118a88cd412fee0c2d0b20c4cce1ed32769f259ac` |
| Accepted cap/asymptotic definitions | `d508cdd3aa82c5204de3626aa24f70199214aa07dc302f3513f246b3eae8fee7` |
| R042 explicit threshold | `a73bae47dbef7b56b95dcc28bf5dc368828265e91c0de49d03a3dcf518e8c0c5` |
| Old R042 comparison/domain chain | `e9dc339e64164f4afb51e20ee81dd03d38f49132353469849de906ec84844c6a` |

Only the independent clone's `cube-dim-7-bounds/audit/round045/comparison/` directory was written for this audit. The parent retains the full-operator proof and final integration responsibility.
