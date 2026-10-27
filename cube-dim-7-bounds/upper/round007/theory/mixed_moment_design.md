# R007 unexecuted matrix-moment candidate for a separate R008 pilot

This is Sol derivation/design only,2026-10-27. No mixed integral, moment scalar, event probability or new I budget is computed here. R007's released numerical activity is only the separately owned matrix-slack diagnostic. Any pilot below requires a new main-owner registration/release in R008; it is not an extension of that diagnostic.

## 1. Target the actual16-column pair matrix moment

Let A,B be two independent raw8-direction block covariance matrices, T_A=det(A)^(1/7), T_B=det(B)^(1/7), L=T_A+T_B, R=det(A+B)^(1/7), H=R/L>=1. Put s=1/8 and q=E L^(-s), already rigorously certified in R004. Under the probability tilt

    dnu=L^(-s)dP/q,

one has the exact identity

    J_pairmat:=E det(A+B)^(-1/56)=q E_nu H^(-1/8).                  (1)

Four independent matrix pairs then give M<=product_j J_pairmat,j/2. A common upper bound J_* for every pair type gives M<=J_*^4/2, without assuming identical pair-matrix distributions. The aim is a material strict upper bound for(1), not another scalar-root AM–GM approximation.

The combined16-column frame is tight with Gram16I. Writing c_J=det(V_J)^2 for7-column subsets, Cauchy–Binet gives

    det(A+B)=sum_(|J|=7)c_J product_(i in J)S_i
      =4^(-7)P9(Y)/product_(i=1)^16Y_i,
    P9(Y)=sum_(|J|=7)c_J product_(i not in J)Y_i.                  (2)

P9 is a positive-coefficient homogeneous degree9 polynomial, at most binom(16,7)=11440 terms. No polynomial or minor has been enumerated in this theory task.

A crucial exact coefficient sum avoids an immediate11440-term integration campaign. If j of the7 selected columns belong to A,

    C_j:=sum_(|J cap A|=j)c_J=8^7 binom(7,j).                      (3)

Indeed these are the coefficients of det(t V_AV_A^T+V_BV_B^T)=8^7(t+1)^7. This uses both8-column frames, not a16-column weighted AM–GM bound. The old16-frame all-AM–GM shortcut returns M_Gamma and remains rejected.

## 2. Exact exchangeable mixed Gamma transform

For a selected set of j coefficients in one block and e>=0, define

    F_j(z;e)=E[(product_(selected j)S_i^e) T^(-z)].

The symmetric8-term determinant identity makes this depend only on j, not which columns were selected. With Y_i=Rp_i, R~Gamma4 and p~Dirichlet(1/2)^8,

    F_j(z;e)=2^(-4z/7)4^(-je)
       ·Gamma(4+z-je)/Gamma(4+8z/7-je)
       ·Gamma(1/2+z/7-e)^j Gamma(1/2+z/7)^(8-j)/Gamma(1/2)^8.      (4)

The derivation is direct: detSigma_block=16R^(-7)/product p, while product S_i^e=4^(-je)R^(-je)product_selected p_i^(-e). The radial power is z-je. Conditions include1/2+Re z/7-e>0 for selected entries and4+Re z-je>0; the unselected parameters must also be positive. A common sufficient domain for j0,...,7 is

    Re z>7e-7/2.                                                  (5)

For j0 this is stronger than necessary, harmless. On this domain F_j is analytic, positive on the real axis, log-convex there and |F_j(v+it;e)|<=F_j(v;e), by its positive weighted-moment representation. Gamma denominators are reciprocal-Gamma factors; they introduce no poles or zero-denominator issue in the stated domain.

For b>0 and a=b/2 define the mixed scalar expectation

    I_j(e,b)=E[(product of7 selected S_i^e)L^(-b)],

with j selected in A and7-j in B. Absolute Mellin/Fourier inversion gives

    I_j(e,b)=1/[pi Gamma(b)] integral_0^infinity Re K_j(t;e,b)dt,
    K_j=Gamma(a+it)Gamma(a-it)
                       F_j(a+it;e)F_(7-j)(a-it;e).                (6)

**The kernel is not generally a modulus square.** Its real part may have either sign and its imaginary part need not vanish at a positive t. The full-line integrand obeys K_j(-t)=conjugate(K_j(t)), which justifies the half-line real part. Only I_j itself is positive. Every signed panel coefficient and signed panel integral must be retained. Pair-exchange symmetry gives I_j=I_(7-j), so only j0,1,2,3 are needed.

## 3. Strong lower bound for E_nu H from four1D integrals

Set e=1/7,b=9/8. For a7-subset J put x_J(S)=c_J^(1/7)product_(i in J)S_i^(1/7). Its l7 norm is R=det(A+B)^(1/7). Minkowski's integral inequality for the positive weight L^(-9/8) gives

    E[R L^(-9/8)]
       >=[sum_J c_J I_j(1/7,9/8)^7]^(1/7)
        =8[sum_(j=0)^7binom(7,j)I_j^7]^(1/7).                    (7)

This is also the optimised fixed-weight root-Jensen inequality: choose constant weights alpha_J proportional to c_J I_j^7, independent of the random S. Thus it is stronger than the natural linear CB-Jensen sum. No data-dependent random weights are substituted inside an expectation.

Consequently, with I_j=I_j(1/7,9/8),

    E_nu H>=L_H:=8[2(I0^7+7I1^7+21I2^7+35I3^7)]^(1/7)/q.        (8)

All lower endpoint operations are monotone: nonnegative complete I_j lower endpoints and the complete q upper endpoint give a safe L_H lower. Define mu_L=L_H_lower-1. The first gate is mu_L>0. If this gate fails, it rejects this particular lower surrogate, not the actual matrix gain or the operator.

The R007 sample hierarchy suggests material matrix gain, but it does not prove(8) exceeds1. That must be checked by the separately released finite pilot.

## 4. Exact seventh moment of H and a finite second-moment bound

The same coefficient identity(3), now with e=1,b=7+s=57/8, gives an **exact** moment

    K7:=E_nu H^7
       =(8^7/q)sum_(j=0)^7binom(7,j)I_j(1,57/8)
       =(2·8^7/q)[I0+7I1+21I2+35I3].                             (9)

All these weighted scalar moments are finite. Their Mellin centre is a=57/16, while the common mixed-moment domain is Re z>7/2. The distance to that pole boundary is only1/16: using a Cauchy radius3a/4 would be **invalid**, despite the Gamma kernel itself having a much wider strip. This is a genuine mixed-transform constraint.

For X=H-1>=0, E_nu X=mu, E_nu X²=nu. Lyapunov and the positive binomial expansion give two strict population upper formulas,

    nu<=K7^(2/7)-1-2mu,
    nu<=(K7-1-7mu)/21.                                           (10)

Therefore an upper K7 and lower mu produce safe upper endpoints for both expressions, and their minimum is a safe nu_U. It must be positive and at least mu_L²; a failure of these enclosure checks is unresolved numerical certification, not a license to clip.

For g(X)=(1+X)^(-1/8), g'''<0. The exact two-moment Hermite/Radau argument from R005 gives

    E_nu H^(-1/8)<=1-delta,
    delta(mu,nu)=mu²/nu·[1-(1+nu/mu)^(-1/8)].                     (11)

Delta increases with mu and decreases with nu. Thus mu_L>0,nu_U>=mu_L² give a safe delta_L by downward evaluation. Require0<delta_L<1. Equations(1),(11) then prove the candidate strict pair-matrix budget

    J_pairmat<=q_upper(1-delta_L),
    M<=q_upper^4(1-delta_L)^4/2,
    I<=M_upper^(7/6)/(105·2^(25/6)pi^(13/12)).                     (12)

The common bounds(8)-(10) use only exchangeability and tight frames, so apply to every matrix pair, even if their full laws differ. This candidate retains a matrix-root-ratio moment, rather than reverting to full16-frame AM–GM. Its quality is unknown until the pilot.

## 5. A cheaper alternative upper moment, if separately useful

For any nonnegative exponents e_i and block determinant exponent B,

    E[det(Sigma_block)^(-B)product_i S_i^(e_i)]
      =2^(-4B)4^(-sum e_i)
        Gamma(4+7B-sum e_i)/Gamma(4+8B-sum e_i)
                      [product_i Gamma(1/2+B-e_i)]/Gamma(1/2)^8.  (13)

For E_nu H², use(subadditivity at power2/7)

    det(A+B)^(2/7)<=sum_J c_J^(2/7)product_(i in J)S_i^(2/7)

and L^(-17/8)<=2^(-17/8)(T_AT_B)^(-17/16). This yields

    E_nu H²<=2^(-17/8)/q
       ·sum_(j=0)^7 C_j^(2/7;sum)
             F_j(17/16;2/7)F_(7-j)(17/16;2/7),                    (14)
    C_j^(2/7;sum)=sum_(|J cap A|=j)c_J^(2/7).

All Gamma parameters are positive. Unlike(3), the powered coefficient sums require an exact minor/histogram enumeration, up to11440 minors per distinct paired frame type. They are algebraic numbers, not powers of C_j. Replacing sum c_J^(2/7) by(sum c_J)^(2/7) would reverse the required bound and is invalid.

Equation(14) is an optional finite Gamma-only second-moment upper bound. Its subadditivity loss may be large. It is not automatically preferable to(9)-(10), and no such enumeration or evaluation is performed here. One may use the smaller of independently certified valid upper bounds under a future fixed registration.

## 6. Signed mixed-integral certification: domain, envelope and full tail

For(6), write H_j(v)=Gamma(v)F_j(v;e). On the common real domain this is positive log-convex. For a complex frequency t+iy, the two Mellin arguments have real parts a-y and a+y. A valid radius rho must satisfy

    a-rho>max(0,7e-7/2).                                         (15)

Log-convexity of H_j(a-y)H_(7-j)(a+y) in y gives the disk envelope

    Cstar_j=max[H_j(a-rho)H_(7-j)(a+rho),
                H_j(a+rho)H_(7-j)(a-rho)].                        (16)

It is the maximum of **two crossed endpoint products**, not a symmetric single product unless j=7-j. There is no assumption of monotone Gamma values.

For T>0,N>=1, Gamma recurrence and the positive weighted-moment bound imply

    |normalised tail beyond T|
      <=F_j(a;e)F_(7-j)(a;e)Gamma(a+N)^2/
                         [pi Gamma(b)(2N-1)T^(2N-1)].             (17)

The sign of this tail is unknown. Subtract its absolute upper bound from the lower endpoint and add it to the upper endpoint. This differs from R004's nonnegative-tail assembly. An I_j lower endpoint can be intersected with0 only because the **original expectation** is known positive; do not replace the integrand by its absolute value when seeking that lower bound.

On panels of halfwidth h<rho, a degree n polynomial has raw aggregate error

    T Cstar_j(h/rho)^(n+1)/(1-h/rho).                              (18)

Integrate every even **real part** coefficient, including negative terms. A panel integral need not be positive. Individual imaginary coefficients need not contain0. A correct series implementation multiplies Gamma(z)F_j(z;e) by the coefficientwise conjugate of Gamma(z)F_(7-j)(z;e), not by its own conjugate. Direct point/series overlaps and full-line conjugacy are useful implementation diagnostics, not a proof that this is a modulus square.

## 7. Proposed finite R008 pilot and materiality gates — not released here

StageI tests only(8): four mixed integrals j0..3 with e=1/7,b=9/8,a=9/16. Suggested fixed parameters are T32,Ntail16,512-bit balls,degree128/order129,rho27/64,halfwidth1/4,64 panels per integral (256 total), eta16/27. All Cauchy-domain conditions(15) hold. Use a4-panel pilot per integral and reuse it. All signed terms and symmetric absolute tails are retained. This is a small original1D Gamma problem with no new shape sample.

After complete StageI, require mu_L>0. Also enclose the explicit surrogate L_H from above. Even under the optimistically smallest feasible variance nu=mu², a certificate based on this surrogate cannot do better than

    M_pair/(1+mu_candidate)^(1/2).

Thus its optimistic lower budget using the surrogate's upper endpoint can be compared with the current R006 expression or a declared material target such as31/500. This is a **limit of this certificate surrogate**, not a lower bound on actual M. If it cannot reach the declared target even optimistically, stop and record that specific design failure before StageII.

Only a new explicit release after these gates can start StageII for(9): four mixed integrals with e1,b57/8,a57/16, suggested T64,Ntail32,512-bit balls,degree192/order193,rho3/64,halfwidth1/32,1024 panels per integral (4096 total), eta2/3. The narrow radius is dictated by the Re z>7/2 mixed pole, not by a alone. Error and full signed tail use(16)-(18) with the actual parameters. No result is assumed from these conservative counts.

Then assemble K7_U,nu_U,delta_L and(12), compare with the complete R006 **budget-expression lower** endpoint, and separately test the declared material target. If the correction is non-material or comparison is unresolved, report failure of this finite surrogate; no automatic additional moments, epsilon search or63D partition is authorised.

All these are unexecuted candidate parameters. R007 supplies no pilot result, strict matrix moment endpoint or new I bound.
