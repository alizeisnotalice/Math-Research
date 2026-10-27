# Matrix-gap subtraction: exact expectation interface and a non-material event test

This is unexecuted Sol derivation, not a numerical budget. It answers precisely when a gap may be subtracted from R006's distribution-moment bound.

Let raw block covariances be Sigma_g, T_g=det(Sigma_g)^(1/7), X=sum_(g=1)^8T_g. Put

    Fscalar=X^(-1/2), Factual=det(sum_gSigma_g)^(-1/14).

Determinant Minkowski gives0<=Factual<=Fscalar everywhere. R006 proves **E Fscalar<=M6**. It does not claim M6 is a per-row majorant. If a deterministic event argument proves E(Fscalar-Factual)>=eta, then

    actual M=E Factual
      =E Fscalar-E(Fscalar-Factual)<=M6-eta.                       (1)

This is the valid subtraction interface. A Monte Carlo mean gap or the unrelated difference Z2-Zpairmat cannot simply be subtracted from M6.

## 1. An explicit finite-probability matrix-pair event

Choose one pair of8-column frames. For any column v in its first frame, tightness of the other frame gives sum_w(v·w)^2=56. Dot products are odd integers. If every |v·w|>=3, the sum would be at least72. Thus each v has at least one w with |v·w|=1. There are at least8 such cross pairs and at most64; an exact64-dot-product check can count N_close if ever needed.

Fix0<ell<u<Llo<Lhi. On an event E_(v,w), put the two selected coefficients S_v,S_w in[Llo,Lhi] and the other14 pair coefficients in[ell,u]. Since Llo>u, events for different selected(v,w) are disjoint. Their total probability is

    p_pair=N_close p_base^14 p_spike²,
    p_base=erf(1/(2sqrt(ell)))-erf(1/(2sqrt(u))),
    p_spike=erf(1/(2sqrt(Llo)))-erf(1/(2sqrt(Lhi))).                (2)

No coefficient independence beyond the original iid S_i is invented.

On E_(v,w), with A,B the two matrices,

    det(A)^(1/7),det(B)^(1/7)<=U0,
    U0=[(8u)^6(u+7Lhi)]^(1/7),
    det(A+B)^(1/7)>=V0,
    V0=[(16ell)^7(1+7t/(8ell)+3t²/(16ell²))]^(1/7),
    t=Llo-ell.                                                    (3)

The upper follows from A<=8uI+(Lhi-u)vv^T. The lower follows from A+B>=16ell I+t(vv^T+ww^T), and the rank2 determinant with |v·w|=1. Hence the root gap G0=V0-2U0 is positive whenever the explicit scalar inequality V0^7>128U0^7 holds.

A guaranteed rational family is u=kappa ell,Llo=tau ell,Lhi=2Llo, with kappa>=1,tau=32kappa^6. The positivity condition reduces to

    kappa^6(kappa+14tau)
       <8+7(tau-1)+(3/2)(tau-1)^2.

Its left side is at most449kappa^12, while the quadratic term alone is at least768kappa^12 since tau>=32. Thus the gap is strictly positive without a numerical search.

## 2. The other six blocks can be controlled by one exact Gamma moment

Let mu_T=E T=m(-1), finite by the block Gamma formula (its domain extends beyond−1). Explicitly

    mu_T=2^(4/7)Gamma(3)/Gamma(20/7)
                         [Gamma(5/14)/Gamma(1/2)]^8.

For K=12mu_T, Markov gives P(sum_other6 T_g<=K)>=1/2. This event is independent of the selected pair. On the combined event, X<=Xmax=2U0+K and determinant Minkowski after merging the selected matrix pair gives

    det(sum_gSigma_g)^(1/7)>=X+G0.

The function x^(-1/2)-(x+G0)^(-1/2) decreases in x, so

    E(Fscalar-Factual)>=eta,
    eta=(p_pair/2)[Xmax^(-1/2)-(Xmax+G0)^(-1/2)]>0.                (4)

Equations(1),(4) give a genuinely strict matrix improvement of the R006 scalar budget. It integrates the original **raw** S law; do not multiply eta by the global radial factor D again. No63-dimensional tiny box is involved.

## 3. Why these coarse events should not be the material R008 route

Consider the four concrete scale bands(ell,u)=(1/4,1),(1/8,1),(1/8,2),(1/16,1), using the guaranteed tau above. They all have u<=2 and Llo>=2^15. They are finite valid strict certificates, but their maximum possible contribution is already analytically too small to justify a new main campaign.

For u<=2, p_base<2/3. One proof is erf(x)>=2(x-x³/3)/sqrt(pi) at x=1/(2sqrt2); with pi<4 this exceeds1/3, so the probability above that lower Y threshold is below2/3. Also p_spike<=1/sqrt(pi Llo). Therefore

    p_pair<=(2/3)^14/(512pi).

The block formula gives T>=16/R_g with R_g~Gamma4, so mu_T>=16/3 and K>=64. The bracket in(4) is at most1/8. Thus eta<=p_pair/16. Meanwhile every valid M6 exceeds actual M>=1/(32sqrt(pi)), by Hadamard's determinant bound and the stable law of sum64 S_i. Consequently

    eta/M6<=(2/3)^14/(256sqrt(pi))<1/12800<10^(-4).                (5)

The last rational inequality follows already from(2/3)^10<1/50 and pi>1. Thus even an optimistic exact evaluation of these four events cannot give a0.01% relative M improvement; their guaranteed contribution will usually be smaller. This excludes the **specified event certificates as material designs**, not the matrix gain itself. The I relative change is at most7eta/(6M6), also microscopic.

Do not spend a new strict-budget round evaluating these events merely because(4) is strictly positive. Their role is to validate the expectation-subtraction logic and explain why a stronger mixed-moment/conditional certificate is needed. The separate mixed-moment design is the prospective material candidate; none of it is executed here.
