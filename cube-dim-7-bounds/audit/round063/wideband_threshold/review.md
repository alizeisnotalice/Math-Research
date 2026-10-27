# R063 independent explicit wide-band threshold audit

All logarithms are natural, as in the accepted operator framework.

Conclusion: **FULL ACCEPT**, for the new R063 domain

    every integer dimension n with log n>=2^33.

On that domain the accepted R059 finite-history mechanism gives

    Cn<=Phi59=F(A_n(t),B_n(t)-(1637/12500)n)
       <n log n-n log log n-n/10.

This is a new explicit threshold established in R063. It does not retroactively change R059's original eventual theorem, computation receipts, or source identity. No minimal-threshold claim is accepted.

Audited frozen source:

    /Users/zhengzhihao/Desktop/T/cube-r030-parameters/cube-dim-7-bounds/upper/round063/wideband_threshold/THRESHOLD_PROOF.md
    SHA256 70fc76d864f233900baadd7f177b7532b654355c3ce3b515a822fe6b65a19133.

The complete source and its verifier were read, but the author verifier was not executed or imported. The actual independent check uses a different, more conservative lower bound on n and different Taylor orders. Only this private clone's audit directory is written; no author file, old grid, accepted verifier, production batch, shared MATLAB or Git operation was used.

## 1. Accepted fixed inputs and scope

The R059 source is final `fe95f369...86bd8f`, its complete audit is `c8aff39b...a65a33d4`, and its actual local lemma is final `97430af9...a243eec`. Their identities were checked. Their already accepted fixed logarithm/moment prices are quoted, not rerun: mL=159477/10^6, Qbar=26979/250000, mean rate<293/20000, variance rate<1/900, gamma=1637/12500, and full scalar slack 39/1212500. They are algebraic fixed-tuple facts once their geometric hypotheses hold; they do not contain an additional unknown n threshold.

All remaining eventual hypotheses of that mechanism are explicitly covered below. The tuple, source/target intervals, weights and original operator are unchanged. The conclusion remains Lebesgue, centered axis-parallel cubes with one common radius in all coordinates, arbitrary L1 input and every positive radius. It is not a bounded-scale or finite-grid operator assertion.

## 2. Independent compression and why it covers the infinite tail

Instead of the author's N(L)=L^8/8!, use the weaker bound

    n=e^L>=L^8/8!>nbar(L):=L^8/2^16,

because 8!<2^16. At L0=2^33, nbar=2^248 exactly. The weaker bound still passes every required probability and geometric guard. This is a new finite certificate, not a replay of the author's 61 checks.

The independent energy-payment bounds also use different orders:

* source cap: Taylor degree 10, giving a positive constant times L^-9;
* tau: Taylor degree 6, giving a positive constant times L^-9/2;
* complete scalar remainder: Taylor degree 4, giving L^-4 and L^-7/2, plus nbar^-1=constant*L^-8.

The global inequality exp(y)>=y^m/m! is valid for every y>0. Combined with the indicated powers, base checks extend to every L>=L0. The upper sqrt(L0)<2^17 is used only to rationally evaluate the decreasing tail functions at L0; it is not falsely imposed on arbitrarily large L.

## 3. Geometry, full source band and active Rbase

The fixed rational domains are checked: 0<d<1/4<1/3, 1/2<qmin<qmax<1, and 1<Rbase<rlo<rhi<5/4. Also x=97L/100 is between L/2 and L, increases, and exceeds one. Thus t=n/x>=nbar/L>16 and t<n for the whole tail. Both target/source annuli are within ||z||<t.

The independent base conditions include

    h<min(inset/2,d/8,(Rbase-1)/2),
    1/(x*nbar^2)<h/2,
    qmin*nbar/L>4rhi,
    nbar>=200000, nbar*d>=1.

Here h decreases like L^-1, the second bound divided by h decreases like L^-16, and nbar/L increases like L^7. Therefore these are full-tail bounds. They give:

* zeta*t=1/(x*n^2)<h/2<d and <1, hence a0=(1-zeta*t)q>0 and kappa>1;
* beta<2, g=d-h+zeta*t*q in (0,d), and g<R/2;
* the whole-band R=q+d-h>1, keeping G/b well defined even for inactive sources;
* a source with any matched target has R>=Rbase+inset-2h>Rbase;
* each complete matching interval [r-d-h,r-d+h] stays inside the source band and has width exactly 2h;
* q_b=r-d+h makes R=r exactly, and the full baseline-to-match interval meets every guard of the accepted d<=1/3 J/AC lemma;
* source rho>=qmin*nbar/L>4rhi>=4r, while target r<2. Thus A/P are disjoint.

The stronger Rbase bound is not applied to inactive sources. Their columns vanish, as required by the frozen R059 kernel definition. No old narrow-band geometry has been substituted.

## 4. All FIVE target failure probabilities are quantified

On each fixed full target face/radius, the baseline pinned coordinate is deterministic and the other n-1 are independent uniforms. Substituting n>=nbar into the accepted baseline mean/variance gives

    Mlog<=6d^3 L^2/(qmin^2*nbar)+dL/(qmin*nbar),
    Vlog<=12d^3 L^2/(qmin^2*nbar).

The new base check verifies Mlog<e/2 and Vlog<(e/2)^2*epsA/8 for e=2^-20. Their terms decrease as L^-6,L^-7 and L^-6, respectively. Chebyshev therefore gives the baseline log-product bad probability below epsA/8 uniformly on every full target face/radius.

Because d<1/4, Kw<=4nd is not automatic and was not discarded. Its exact law is 1+Bin(n-1,2d/r), with mean<=1+2nd and variance<=n/4. Once nd>=1, the threshold gap is at least nd; hence the bad probability is at most 1/(4nd^2)<=1/(4*nbar*d^2)<epsA/8. This bound decreases like L^-8.

At either count endpoint, the mean is n+(n-1)a0/r and variance<=n/4. The two failures at tolerance e*n/4 are each at most 4/(e^2*nbar)<epsA/8. Endpoint calibration retains the 1/n drift: the independently verified full-tail conditions are `2h+1/nbar<3e/4` and `2hd<3e/4`. This gives the full-q normalized count in [1-e,1+e]. These drifts decrease, not merely their separate point samples.

For the target D cap, the check `1/rhi-L/nbar>1/2` makes every positive-face coordinate fraction `1+1/(2t)-1/(2r)` strictly below 3/4. The negative pinned face contributes zero. Since n>=65, per-face complement probability is at most (3/4)^64<epsA/8. No averaged factor 1/2 is used.

Thus baseline, Kw, count-minus, count-plus and target cap each cost at most epsA/8. Their union costs at most 5epsA/8<epsA; no independence of these events is needed. The resulting target mass is at least (1-epsA)n*log(rhi/rlo)>mL*n by the already accepted exact fixed logarithm bound.

## 5. Product and individual-J guards are separate

On Kw<=4nd, the accepted local a.e. derivative/absolute-continuity bound summed over the n coordinates is at most 640Ld. The full matching interval length is at most 2h, so the log shift is 1280Ldh=1280h_c. Baseline tolerance e therefore supplies the accepted fixed product bracket P0=1-e-1280h_c.

Separately, the pointwise local coordinate estimate is

    |J-1|<=5d/(8rho)<=5dL/(8qmin*nbar)<e.

Its upper bound decreases like L^-7 and is uniform in coordinate value and matching q. It is not a probabilistic event and needs no union over n coordinates. It is not raised to the n-th power in place of the target product condition. These deterministic inputs and the endpoint counts give the exact Ulo/Uhi full signed-source-face row estimates for all L in the claimed domain.

## 6. Finite inverse retention, including its tighter threshold

The active-source geometry now satisfies every domain needed by the accepted fixed inverse moments. Each captured signed source face has

    E Phi<=1/9+(293/20000)x,
    Var Phi<=x/900.

The pin is deterministic, free coordinates are independent only within that full captured source face, and no cap-conditioned source law is used. Since chi-293/20000=1/20000 and x>=40000/9, the cutoff is at least x/40000 above each face's mean. Chebyshev gives

    tail<=40000^2/(900x).

The new exact base calculation yields

    40000^2/(900x0)=390625/1830813696<1/4096=epsRet.

This is among the tighter threshold checks and was independently evaluated; it was not inferred from a generic sufficiently-large statement. It decreases like L^-1. Face probabilities may therefore be mixed while retaining the same epsRet bound, without assigning the within-face variance to the mixture. The full soft-retention row bound holds at every L>=L0.

## 7. Wide-source cap: retain n-1 and the correct Csrc

The source-band mass coefficient is Csrc=(qmax-qmin)/qmin, not d; for this tuple Csrc>d. Its accepted cap upper bound divided by k is

    (Csrc/h)*exp[-((1-1/n)eta-chi)x].

The fixed tuple has 0<eta<1 and eta-chi>1/100000. The new bound n>=nbar>=200000 gives

    (1-1/n)eta-chi
       >=eta-chi-eta/nbar>1/200000.

This preserves the finite n-1 correction. Put cS=97/20000000. Then the ratio is at most `(Csrc/h)exp(-cS L)`. Using Taylor degree 10 gives a full-tail upper bound proportional to L^-9. At L0, Csrc<1 and d<1/4 give Csrc/h0<2^79, while cS*L0>2^15. Consequently the independent bound is

    cap_ratio<10!/2^71<epsS.

This check is stronger than needed but uses a different proof compression from the author's Taylor-eight certificate. The original source deletion is paid by its full unweighted mass because b<=1. No coordinate independence after Dcap is reused.

## 8. Positivity of the selected root and payment of tau

The original global integral estimate I_n(t)<=2t/(n-1), beta<2 and n/(n-1)<=2 imply

    B_n(t)>=n log t-2n/x>=T-4n/L.

For L>=16, log L<L/4: it holds at 16 using log2<1, and the difference L/4-logL is increasing. The independent new base inequality `L/4-4/L-1>0` then gives B_n-gamma*n>nL/2>0 throughout the domain, since gamma<1. On the upper side x>4 and n>=2 ensure B_n<=1+n logt<=nL; A_n<=exp x is exact. Therefore the specified real positive tau is well defined.

Its normalized row payment is at most `(sqrt L/h)exp(-nu L)` with the exact positive nu=741/10^6. Taylor degree 6 gives a full-tail upper bound proportional to L^-9/2. At L0, sqrt(L0)<2^17, d<1/4 and h_c=2^-48 imply sqrt(L0)/h0<2^96, while nu*L0>2^22. Thus

    tau/k<6!/2^36<epsT.

The source cap and tau use separate positive exponents and separate epsilon payments. Both are now quantified; neither follows merely from a positive scalar model objective.

## 9. Every actual finite history and the original operator

For every algorithmic dilation mathfrak a>=8n^4, t<n and n>1 give mathfrak a>2t and c_cap>1. Beta<2 gives the actual displacement guard

    alpha*beta-kappa-2/(mathfrak a+1)
      >n^-3-4/(8n^4)>0.

The original corner-envelope inclusion is algebraic for m>=1 and beta>1, so its exact equality lambda_a=alpha^-n*lambda on A/P has no further unquantified n guard. Actual target items meeting A have radii less than alpha+rhi<3; actual source items have radii at least rho>4. Thus all included sources precede the same target in the original radius-priority history and really capture its center. Since b<=1, integrating included source pieces pays no more than those original item weights, giving the same true eligibility bound tau.

The accepted energy identity uses the same nonnegative missing measure lambda_a-mu_a on disjoint A/P. Ambient energy scales by alpha^-2n and degrees/target mass by alpha^-n. Since alpha^-n>=1, the quantified row excess, target mass and degree budgets give at least Gamma*n>gamma*n missing envelope mass uniformly for every such dilation and every original finite family. No selector or history limit is required.

The old cover factor, captured tau and terminal self term remain in the finite witness inequality. Only those uniform scalar inequalities are passed to mathfrak a->infinity, giving the full root cost Phi59. General L1 and all positive common radii are preserved by the same compact strict-level-set finite cover, actual witness dilates and terminal-radius normalization. Inner regularity yields the full level set. No hidden compact-support, radius-grid, upper/lower scale cutoff or alternative rectangular operator enters.

## 10. Complete scalar remainder, not only the deficit

The accepted exact price is gamma-sigma/(1-sigma)=1/10+39/1212500. The whole positive remainder divided by n is bounded by

    1/n+exp(-sigma L)+2sqrt(L)exp(-sigma L/2).

Our independent Taylor-four upper envelope is a positive sum of powers L^-8,L^-4,L^-7/2. At L0, use 1/n<=1/nbar=2^-248, sigma*L0>2^27, sigma*L0/2>2^26, and sqrt(L0)<2^17. The remainder is strictly below

    2^-248+24/2^108+24/2^86<39/1212500.

All terms decrease on the complete tail. Thus Phi59<T-n/10 for every claimed dimension, with all original cover costs paid.

## 11. Comparison of the two complete root costs on their actual common domain

The R056 accepted coefficient gamma56 was read from its actual accepted receipt, not recomputed. Exact data checks give gamma56<1/8. For its different x56=93L/100, the original integral lower bound and nonnegative core give

    B56>=T-log(93/100)*n-4n/L
        >=T+(7/100)n-4n/L.

Therefore Phi56>=B56-gamma56*n. At L>=2^20, the bound

    gamma56-7/100+4/L
      <1/8-7/100+4/2^20<1/10

shows Phi56>T-n/10. The positive R056 B-buffer 7/100 is essential and retained. Since R063's L>=2^33 lies inside R056's actual L>=2^20 domain, this yields Phi59<T-n/10<Phi56 on the true common domain. It compares these two complete internal functions, not only their displayed negative terms and not all possible strategies or literature records.

## 12. Actual verification and identities

`check_compressed_threshold.py` was actually run once with standard-library Fraction/integer arithmetic. Twenty-five compressed new guard groups and three complete-cost comparison groups passed. It uses nbar=L^8/2^16 and Taylor orders 10/6/4, not the author's 61-check code, and imports no author module. Accepted historical JSON is read only for fixed data compatibility and gamma56; no log series, old budget verifier, grid, production batch or MATLAB computation is replayed.

The source-lock addendum records the final audit/source/receipt identities and the accepted exact gamma56 identity. All original R059 and R056 histories remain unchanged. This independent explicit threshold is attached to R063 only.
