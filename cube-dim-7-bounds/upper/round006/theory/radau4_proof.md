# R006: four-moment Gauss–Radau refinement with strict R005 comparison

2026-10-27. Sol theory only. This owner runs no numerical script, quadrature, MC or MATLAB. The original n7 body64/equal weights/common Cauchy time and the8 disjoint Hadamard blocks are unchanged. This file does not authorise a numerical start; only the main owner may release the separately registered finite batch. No g/C7 claim follows from an I budget alone.

## 1. The same target under a finer log-ratio variable

For T=det(Sigma_block)^(1/7), define m(z)=E T^(-z) for Re z>0 by the already proved formula

    m(z)=2^(-4z/7)Gamma(4+z)/Gamma(4+8z/7)
                    [Gamma(1/2+z/7)/Gamma(1/2)]^8.

Let U=T1+T2, q(r)=E U^(-r), q(0)=1 exactly. The complete R004 pair moment q(1/8) and the three R005 pair moments q(1/16),q(3/16),q(1/4) are read-only reusable inputs. q0 is probability normalisation, not a Mellin integral at0.

Fix

    a=1/8, s=1/4, h=1/32,
    dP_a(u)=u^(-a)dP_U(u)/q(a).

For independent U1,U2 under this tilted law put Delta=log U1-log U2 and

    W=cosh(h Delta)-1>=0,
    g(w)=[T16(1+w)]^(-1/4).

Because1/(2h)=16, cosh(Delta/2)=T16(cosh(h Delta)) exactly. Therefore the identical four-block target satisfies

    q_four(1/4)=E(T1+...+T4)^(-1/4)
       =2^(-1/4)q(1/8)^2 E_a g(W).                                (1)

This changes only the auxiliary variable; the target, tilt and exponent do not change.

The U law is continuous with full support(0,infinity): conditional on an internal block angle, T is a positive constant divided by a Gamma4 radius with full support. Convolution, positive tilting, log differences and cosh then imply W has continuous full support on(0,infinity), with no atom at0. This ensures strict positivity of the moment matrices and quadrature weights below.

## 2. Four exact moments and only four new real pair integrals

Put R_l=q(a-lh)q(a+lh)/q(a)^2, l=0,...,4, with R0=1. They equal E_a cosh(lh Delta). The moments m_j=E_a W^j, j=0,...,4, are

    m0=1,
    m1=R1-1,
    m2=R2/2-2R1+3/2,
    m3=R3/4-3R2/2+15R1/4-5/2,
    m4=R4/8-R3+7R2/2-7R1+35/8.                                   (2)

These follow from expanding(cosh-1)^j and the identities for cosh^j. Equivalently E cosh^j(hDelta)=2^(-j)sum_(k=0)^j binom(j,k)R_|j-2k|.

Here

    R1=q(3/32)q(5/32)/q(1/8)^2,
    R2=q(1/16)q(3/16)/q(1/8)^2,
    R3=q(1/32)q(7/32)/q(1/8)^2,
    R4=q(0)q(1/4)/q(1/8)^2.

Thus the **only** new pair exponents are1/32,3/32,5/32,7/32. All endpoints are positive except the known exact q0. Every moment through4 is finite. This construction never evaluates a complex-parameter q.

## 3. Exact positive nodes from the true moments

Define

    D=m1m3-m2²>0,
    alpha=(m2m3-m1m4)/D,
    beta=(m2m4-m3²)/D,
    p2(w)=w²+alpha w+beta.                                         (3)

The equations are

    m3+alpha m2+beta m1=0,
    m4+alpha m3+beta m2=0.

They state that p2 is orthogonal to1,w for the positive weighted measure w dP_a(W). D is strictly positive by strict Cauchy–Schwarz for that measure. All these integrals require only m1,...,m4.

Here is a self-contained root argument. If p2 had at most one sign-changing zero on(0,infinity), multiply it by the product of those zeros' linear factors. That multiplier has degree at most1, while the product has a fixed nonzero sign on(0,infinity), since every remaining positive real zero has even multiplicity and nonpositive/complex factors cannot change sign there. Its integral against w dP_a is then strictly nonzero by full support, contradicting the two orthogonality equations. Consequently p2 has two distinct positive roots

    0<d1<d2.

In particular alpha<0, beta>0 and discriminant alpha²-4beta>0. This proof does not assume a fifth numerical moment or an unverified interval eigenvalue computation.

## 4. Positive weights and exactness through degree4

For nodes d0=0,d1,d2, let L_i be their degree2 Lagrange cardinal polynomials and define lambda_i=E_a L_i(W). The weights may be evaluated explicitly as

    lambda1=(m1d2-m2)/[d1(d2-d1)],
    lambda2=(m2-m1d1)/[d2(d2-d1)],
    lambda0=1-lambda1-lambda2.                                    (4)

For any polynomial P of degree<=4, divide by w p2(w):

    P(w)=w p2(w)ell(w)+r(w), degree ell<=1, degree r<=2.

Orthogonality makes the first term's expectation0 and it vanishes at all three nodes. The definition of the L_i integrates r exactly. Hence

    E_a P(W)=sum_(i=0)^2 lambda_i P(d_i), degree P<=4.              (5)

Apply(5) to L_i²: lambda_i=E_a L_i(W)²>0, because a nonzero polynomial cannot vanish on the full support of W. Also sum_i lambda_i=1. This proves positivity from the exact distribution, rather than inferring it from an approximate interpolation polynomial.

For interval consistency one may also use

    lambda0=det H2/(m2m4-m3²),
    H2=[[1,m1,m2],[m1,m2,m3],[m2,m3,m4]].                           (6)

H2 is positive definite because it is the Gram matrix of1,W,W² for a full-support positive measure; its determinant and the denominator in(6) are strictly positive. Formula(6) is the Christoffel/linear-algebra identity, or follows by inserting(3) in E[p2(W)/beta]. It can sharpen the interval from1-lambda1-lambda2 by intersection; it must not be used to replace a failed positivity check with a midpoint.

## 5. Global Hermite majorant and the certified expectation

The16 zeros of T16 lie strictly between−1 and1. Consequently

    g(w)=2^(-15/4)product_(j=1)^16
               [w+1-cos((2j-1)pi/32)]^(-1/4)

is completely monotone on[0,infinity). Products preserve the alternating derivative signs. In particular g^(5)<0, g(0)=1, g is positive strictly decreasing.

Let P4 be the unique degree<=4 Hermite polynomial with values at0,d1,d2 and first derivatives at d1,d2 matching g. The global remainder for every w>=0 is

    g(w)-P4(w)=g^(5)(xi) w(w-d1)²(w-d2)²/5!<=0.                   (7)

It applies also outside the node interval; xi belongs to the interval spanned by0,d1,d2,w. At the nodes there is equality. By(5),

    E_a g(W)<=B6:=lambda0+lambda1g(d1)+lambda2g(d2),
    0<B6<1.                                                       (8)

This is the exact four-moment Radau expression. Numerical interval arithmetic only has to enclose the **true** moment-defined roots and weights in(3)-(4); it does not have to prove that a polynomial formed from arbitrary interval or midpoint nodes is a global majorant. The latter would be an invalid substitute for(3)-(8).

The resulting original-process budgets are

    M<=q(1/8)^4 B6²/2,
    I<= [q(1/8)^4 B6²/2]^(7/6)/
                           [105·2^(25/6)pi^(13/12)].               (9)

## 6. Strict improvement over the exact R005 Radau expression

Let Y denote the old R005 variable. Since its step was1/16=2h,

    Y=cosh(Delta/16)-1=2W(W+2).

Let f(y)=[T8(1+y)]^(-1/4). Then f(2W(W+2))=g(W) by T8(T2(x))=T16(x). The tilt and target exponent are the same in both rounds.

For the **true** old moments mu_old=E Y,nu_old=E Y², the R005 quadratic Hermite majorant P_old satisfies P_old>=f on all y>=0 and has equality only at y=0 or y=d_old=nu_old/mu_old>0; strict negative third derivative and its remainder prove the strictness elsewhere. The composite

    Q(w)=P_old(2w(w+2))

has degree<=4 and dominates g. On w>=0 it has only two contact points:0 and the unique positive w* solving2w*(w*+2)=d_old. At least one of the two distinct positive new nodes d1,d2 differs from w*. All its quadrature weights are positive. Therefore

    B6=sum_i lambda_i g(d_i)
       <sum_i lambda_i Q(d_i)
        =E Q(W)=E P_old(Y)=B5.                                   (10)

Here B5 is the exact true-moment R005 Radau expression. Thus the new analytical budget strictly improves it. This does not yet provide a numerical size or a new rational endpoint; those require the finite certificate.

For a same-moment numerical comparison, the old moments can be reconstructed without new integration:

    mu_old=4m1+2m2,
    nu_old=16m2+16m3+4m4.                                         (11)

These identities are also a useful consistency check against the R005 reused pair moments. A rounded upper endpoint of B6 must be compared with a certified lower endpoint of the old deterministic **budget expression**, never with a lower endpoint falsely labelled as actual I or M.

## 7. Safe interval assembly contract

All m_j are outward intervals built from the complete reused/new q endpoints, including their analytic tails. Do not clip negative lower m_j to0, repair PSD by moving a midpoint, or call a near-zero discriminant positive without separation.

Required numerical gates:

1. q endpoints positive; moment intervals finite; D=m1m3-m2² has positive lower endpoint. Compute alpha,beta by(3); certify alpha upper<0,beta lower>0 and discriminant lower>0.
2. Use stable root formulas

       d2=(-alpha+sqrt(alpha²-4beta))/2, d1=beta/d2.

   They enclose the true ordered roots. Require d1 lower>0 and d1 upper<d2 lower. All denominators must exclude0. An interval quadratic formula is used as an enclosure of exact moment-defined roots, not as a definition of arbitrary altered nodes.
3. Evaluate(4) outward, optionally intersect lambda0 with(6). Require all three weight lower endpoints positive, sum interval contains1, and numerical moment residuals through degree4 contain0. These residuals are implementation diagnostics; exactness was proved analytically in(5).
4. A particularly safe expectation upper endpoint uses positivity and total mass1:

       delta6_L=lambda1_L[1-g(d1_L)]
                    +lambda2_L[1-g(d2_L)], rounded downward,
       B6_upper=1-delta6_L, rounded upward.

   Since g decreases, each correction is a true lower bound. Require0<delta6_L<1. The strict upper bound1 follows also from lambda1_L+lambda2_L<1. Full direct interval evaluation of(8) can be intersected if it is sharper.
5. With q_a_upper from the complete R004 certificate, form

       M_upper=(q_a_upper)^4 B6_upper²/2,
       I_upper=M_upper^(7/6)/(105·2^(25/6)pi^(13/12)).

   Every whole expression is outward. Certify quantitative improvement against the complete accepted R005 budget-expression lower endpoint (and, when convenient, B6_upper<B5_lower from(11)). Then round I upper upward to a declared rational grid, retaining the finer exact dyadic endpoint.

Evaluate T16 by four Chebyshev doubling steps at x>=1, and use positive factorizations of T_(2n)-1=2(T_n-1)(T_n+1) to avoid cancellation. With R=T16^(1/4), 1-g=(T16-1)/[R(R³+R²+R+1)] is an equivalent stable expression. No negative weight or uncertain root is silently clipped.

## 8. The four new Mellin certificates and parameter-specific errors

For each new r, set a_r=r/2 and reuse the proved positive integrand

    f_r(t)=|Gamma(a_r+it)m(a_r+it)|²,
    q(r)=1/[pi Gamma(r)]integral_0^infinity f_r(t)dt.

Take T=32,512-bit outward acb_series, degree128/order129, tail recurrence N=16. The actual fixed grids are:

|r|a_r|rho_r=3a_r/4|halfwidth h_r|panels|eta=h_r/rho_r|
|---|---|---|---|---:|---|
|1/32|1/64|3/256|1/128|2048|2/3|
|3/32|3/64|9/256|1/64|1024|4/9|
|5/32|5/64|15/256|1/32|512|8/15|
|7/32|7/64|21/256|1/32|512|8/21|

All h_r<=a_r/2 and rho_r<a_r. Exact centres are(2j+1)h_r. Total new panels4096; R004/R005 moments are reused, not recomputed.

The parameter-specific Cauchy bound is

    Mstar_r=Gamma(a_r-rho_r)Gamma(a_r+rho_r)
                          m(a_r-rho_r)m(a_r+rho_r).

Log-convexity of Gamma(v)m(v) proves it for the entire disk as in R004. Here a_r-rho_r=a_r/4 and a_r+rho_r=7a_r/4. The full raw finite-segment Taylor error is

    Eraw_r=32Mstar_r eta^129/(1-eta),                              (12)

not a universal copied1/8 error. The normalized full tail is

    Etail_r=m(a_r)^2Gamma(a_r+16)^2/
                          [pi Gamma(r)·31·32^31].                 (13)

Integrate all even coefficients with2h_r^(2k+1)/(2k+1), retain all signed terms, divide finite bounds bypi Gamma(r), and add(13) only to the q upper bound. First omitted Taylor degree is129. The three smaller eta ratios are used as registered; unused CPU is not grounds to enlarge the mesh.

## 9. Finite acceptance and stopping

The accompanying `registration.md` is a proposed finite batch only, pending explicit main-owner release. First lock all source/config hashes and exact reused moment endpoints. Use4 pilot panels for each new r (16 total), then reuse them if the pilot passes. Upon release,15 single-thread processes compute the remaining4080 fixed panels, with actual resource receipts. No MC, point optimiser or extension follows automatically.

Assembly must verify4096 unique(r,index) keys, per-r parameters/degree/precision, unchanged sources, complete tails, moment/root/weight gates, correction and strict comparison. If separation fails, report the precise gate and stop without increasing degree/precision/cutoff or replacing intervals by midpoints. Full analytic strictness(10) remains true, but an unseparated finite certificate is not a new numerical result.

This theoretical task ends at the frozen proof/registration and explicit handoff. It contains no executed mathematical computation and no new g or C7 certificate.
