#146 general-cost interface: independent theory preflight

**PASS for the general-cost theorem and the specified raw rational candidate construction.** This is not yet a review of a frozen146 driver or an executed H certificate. No new H evaluation, optimizer, old objective/baseline, native kernel or moment integration is run in this review. A bounded exact Fraction check of the saved coefficients, cap sources and target formulas passed.

Use the same original body32 mixture as123/125/128: d=det(Sigma)^(-1/2), M=E[d^(1/6)], Phi_0=1 and Phi_(h,k)=Phi_k(N_h). Write mu_i=E[d Phi_i] and let beta_i be **legal upper caps**, mu_i<=beta_i. For any fixed finite nonnegative coefficient vector a with a0>0, define

    W_a=sum_i a_i Phi_i,
    H(a)=E[W_a^(-1/5)],
    C(a)=sum_i a_i beta_i.

This C is an expectation budget, not the actual weak-type model cost. Since W_a>=a0>0, H is finite; C is finite and positive. Linearity gives E[d W_a]<=C without any independence assumption between d and the nested counts. Hölder with exponents6 and6/5, applied to d^(1/6)=(d W_a)^(1/6)W_a^(-1/6), proves

    M <= E[dW_a]^(1/6) H(a)^(5/6),
    M^6 <= E[dW_a] H(a)^5 <= C(a) H(a)^5.             (1)

Thus an accurate complete H_U yields the valid complete M upper budget C(a)H_U^5. The cost appears to the first power, not the fifth. Nonnegativity and a0>0 must be checked for the actual saved rational a; the optimistic normalization cannot replace this legal-cost calculation.

The construction is homogeneous in coefficients only: for lambda>0, W_(lambda a)=lambda W_a, C(lambda a)=lambda C(a), H(lambda a)=lambda^(-1/5)H(a). Hence C H^5 is invariant. Setting a'=a/C(a) gives legal cap normalization C(a')=1 and H(a')=C(a)^(1/5)H(a), exactly recovering(1). This does not change an auxiliary h, the count law or the physical Cauchy time; it is not a spatial/time rescaling.

For the144 candidate let x be its exact nonnegative2^-44 simplex and let tilde_beta agree with the128 legal cap vector except at flat indexj=32, meaning h-index5(h=1/128),rank2. Define **exactly** a_i=x_i/tilde_beta_i, with tilde_beta_j=bL. Then sum_i a_i tilde_beta_i=1 is merely an optimistic algebraic normalization. For a legal cap b at that index and unchanged legal caps elsewhere,

    C(a)=sum_(i!=j)x_i+x_j b/bL
        =1-x_j+x_j b/bL.                              (2)

The saved144 coefficient is x_j=10680174170037/17592186044416. Its x, not generally its quotient a, is dyadic: the a_i are exact rational quotients. The146 author confirms no subsequent dyadic quantization of a. If coefficients are ever rounded, (2) no longer follows automatically; recompute sum actual_a_i*legal_beta_i and all resulting W tables instead.

The checked current legal b is128 upstream beta_by_h[5][1], also124's old upper and130's old_beta2_interval[1]:

    b=40485654615133438144297492278807668644266365541005757703790265979
      /220855883097298041197912187592864814478435487109452369765200775161577472.

The bL in144 is exactly130's true-moment lower endpoint moment_interval[0], approximately2.2443796428811091e-8. The130 true-moment upper is approximately4.669488925556094e-7, larger than the above legal old cap b≈1.8331254774543401e-7; therefore taking the smaller old cap is valid and stronger among these available certified uppers. The exact source check found0<bL<b<130_upper. **bL is not a legal cost upper cap**: its only role here is defining a. The direct sum of all49 actual rational a_i times the legal128 caps equals(2) exactly and is C≈5.3514489795759976>1. This is a saved-coefficient/cap identity, not a new H or M result. Only the changed beta entry differs between144 and128; beta0 and the other47 band caps, h, Phi and CDF inputs agree exactly.

The125 full-count certificate mechanism continues to apply to H for arbitrary fixed nonnegative a: whole-CDF coupling precedes endpoint-law relaxation, all735471 prefix states retain their true positive probability, and the16 K=17..32 tails use the actual pointwise W-tail lower. The coefficient-scale/root checks must enclose the actual raw W_a. An optimistic beta need not be a moment cap to define that positive W, but it cannot justify E[dW]<=1. An H_L, even with C H_L^5 above a threshold, bounds the fixed cap-budget expression from below, **not actual M^6**. A noncompetitive new upper excludes neither actual M nor the family.

By123 the complete transfer is

    I^5 <= (3/16)^6 C(a) H_U^5 /384^5,
    U_certificate(g)^5=(9/8)^6 C(a)H_U^5/(384^5 g^6).

For a separately certified global floor g, the strict below64 test is C(a)H_U^5<Theta(g), where Theta(g)=(8/9)^6*24576^5*g^6. The checked exact current targets are

    Theta(1/65000)=36893488147419103232/629200415790081024169921875,
    Theta(1/63500)=2361183241434822606848/35005089051913616180419921875.

The complete128 saved M_U^6 is an upper certificate, not an actual-M lower endpoint. Exact source comparison places it above the first target and below the second; that neither establishes global gstar nor makes a failure of a new candidate an actual-M obstruction. Compare146's new C H_U^5 against that full saved128 rational upper, never a displayed decimal or144 float objective. Both floor targets remain conditional until their global coverage is proved. No global g or C6 breakthrough is asserted here.

Checked scientific source hashes:

-144 optimistic_inputs.json:167de796551fea26529e6430ba1f0c074e2f5085f78d4f252fe932badbc5969b
-144 result.json:cb2441dbd0d52159448302acb1d0ad0a5b10d9bc0585d99d2dca6565408bcd84
-128 upstream.json:850eda5b98f73b2f50e744cfe5a6ef7e8104a22dbf7decd6797fe22dd7448f37
-128 candidate_certificate.json:15ee2f6d0c978500e67b7267a657984e802d269b44066bcf2a7597f41ea7d301
-130 results.json:e7d3b1ce50a0fdb834734e69d406110cb7e2364639822b561577904520c0f3cb

Mathematical sources read:123 shape_Minkowski_moment/theorem.md §§1,2,5;125 joint_count_moment/proof.md §§1–3,5;128 proof/upstream and the saved full certificate;144 PROOF and exact saved candidate. Final146 frozen source, fixed precision/work partition and saved integer/root certificate still require their own static/independent gate. The main agent retains final responsibility.
