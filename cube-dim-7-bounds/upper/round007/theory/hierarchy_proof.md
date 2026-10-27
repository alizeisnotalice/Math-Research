# R007: per-shape matrix/scalar slack hierarchy and uniform envelope

2026-10-27. Sol proof only, no numerical execution by this owner. Original7D body64, weight1, common Cauchy time and the same8 disjoint Hadamard blocks are fixed. The numerical diagnostic is not a strict I certificate and proves no g or C7 result.

## 1. Exact radial conditioning

Let independent Y_i~Gamma(1/2,1), S_i=1/(4Y_i), Sigma=sum_iS_i v_iv_i^T. Put R=sum Y_i and P_i=Y_i/R. Then R~Gamma32 is independent of P~Dirichlet(1/2)^64, and

    Sigma=A(P)/(4R), A(P)=sum_i v_iv_i^T/P_i,
    D=2 E sqrt(R)=2Gamma(65/2)/Gamma32.

Every function below is homogeneous of degree−1/2 in the7th-root determinants of its matrix arguments. Thus its raw-S counterpart is multiplied by2sqrt(R), and exact conditioning supplies the common factor D. No fresh R draw is needed to compute any observable.

For block g let B_g=sum_(i in g)v_iv_i^T/P_i, r_g=det(B_g)^(1/7). Pair blocks(1,2),(3,4),(5,6),(7,8), and put C_j=B_(2j−1)+B_(2j), R_j=det(C_j)^(1/7). All these matrices are SPD for every open-simplex P.

Define the six nonnegative observables

    Z0=D det(A)^(-1/14),
    Z1=D (sum_(g=1)^8 r_g)^(-1/2),
    Z2=(D/2) product_(j=1)^4(r_(2j−1)+r_(2j))^(-1/8),
    Z3=(D/sqrt(8))product_(g=1)^8 r_g^(-1/16),
    Zpairmat=(D/2)product_(j=1)^4 det(C_j)^(-1/56),
    ZpairSum=D(sum_(j=1)^4 R_j)^(-1/2).                            (1)

The matrix-pair exponent1/56 is correct: four-group AM–GM gives exponent−1/8 on R_j=det(C_j)^(1/7).

## 2. The diamond, not an incorrect total order

SPD determinant Minkowski gives

    det(A)^(1/7)>=sum_jR_j>=sum_gr_g,
    R_j>=r_(2j−1)+r_(2j).

Ordinary AM–GM applied to the indicated four or eight scalars now proves

    Z0<=Z1<=Z2<=Z3,
    Z0<=ZpairSum<=Z1,
    ZpairSum<=Zpairmat<=Z2.                                       (2)

In particular Z0<=Zpairmat<=Z2. **There is no general order between Zpairmat and Z1.** The first replaces within-pair matrix slack but pays AM–GM across four matrix pairs; the second keeps the exact sum of eight scalar roots. The diagnostic must not assert an invented ordering between them or clip a signed difference to0.

Useful nonnegative paired gaps are

    Z1−Z0, ZpairSum−Z0,
    Z1−ZpairSum, Zpairmat−ZpairSum,
    Z2−Z1, Z2−Zpairmat, Z3−Z2,
    Zpairmat−Z0.

The signed difference Zpairmat−Z1 distinguishes the incomparable branches. All gaps are evaluated on the same P row, not from independent sample means.

## 3. Every observable is bounded by D/64

Let beta_g=sum_(i in g)P_i and p_gi=P_i/beta_g. Then B_g=A_g(p_g)/beta_g, where the7x8 block frame H has HH^T=8I and all squared7-minors equal512².

Cauchy–Binet and the uniform8-minor weights give

    det A_g=8^7 (1/8)sum_(J size7)product_(i in J)p_gi^(-1)
       >=8^7 product_(i=1)^8 p_gi^(-7/8)
       >=8^14.

The first AM–GM uses that each column belongs to7 of the8 minors. The second uses product p_gi<=8^(-8). Thus

    r_g>=64/beta_g.

Since sum beta_g=1 and product beta_g<=8^(-8),

    Z3<=D/sqrt(8)·64^(-1/2)(product beta_g)^(1/16)
       <=D/sqrt(8)·(1/8)·8^(-1/2)=D/64.                          (3)

All other observables are below Z3 by(2). Therefore every six-level observation and every nonnegative gap belongs to[0,D/64]; the signed difference Zpairmat−Z1 belongs to[−D/64,D/64]. This is an exact per-P statement, including arbitrarily small positive P_i. Uniform P attains the common upper value.

## 4. What their expectations mean

Exact radial conditioning gives E_P Z0=actual M. The remaining expectations are successive legitimate upper bounds, but a Monte Carlo mean is not a certified value of any of them:

- E Z1 is the exact eight-scalar-root determinant-Minkowski budget.
- E Z2 is the exact R004 pair-AM–GM budget q(1/8)^4/2.
- E Z3 is the R002 full eight-block AM–GM Gamma expression.
- E ZpairSum is the exact four-matrix-pair determinant-Minkowski sum budget.
- E Zpairmat is the four-matrix-pair AM–GM budget. In raw-S variables it is(1/2)product_j E det(Sigma_(2j−1)+Sigma_(2j))^(-1/56), since pairs are independent. Equal pair distributions require a separately checked row-sign/column-permutation equivalence of the paired16-column frames; they must not be assumed merely from disjointness. The per-P proof does not require them to be identically distributed.

The small numerical change from R005 to R006 alone does not show that the remaining scalar slack E Z2−E Z1 is small: the Radau relaxation may itself still be loose. This diagnostic directly compares that gap with E Z1−E Z0 and with the matrix-pair gains. It is the correct gate before more scalar moments or a matrix-pair integration campaign.

R006's rigorous scalar budget M6 satisfies E Z1<=M6. It is a distribution-moment bound, not a per-P constant majorant. Any strict matrix subtraction must prove an expectation gap, as specified in the separate matrix-candidate derivation; one cannot subtract an observed sample gap from M6.

## 5. Numerical evaluation and what the bound permits

A stable implementation uses log determinants and log-sum-exp. All frames used here are tight: the block, matrix-pair and full matrices have unweighted Grams8I,16I,64I respectively. For positive coefficient weights w_i=1/P_i,

    c min_i w_i I <=sum_iw_iv_iv_i^T<=c max_iw_i I

on each such frame, c=8,16 or64. Hence every matrix condition number is at most max_iP_i/min_iP_i on its columns, and the global64-coordinate ratio bounds all of them. This is a cheap mathematical condition upper bound, not an SVD guess.

The numerical contract flags a row for same-input high precision if the global Y_max/Y_min>2^40, any Cholesky/nonfinite failure occurs, or any asserted hierarchy/range invariant fails beyond the declared arithmetic tolerance. Since Y ratios equal P ratios, no extra solve is needed for that flag. Unflagged condition bounds alone are not a rigorous total floating-error proof; this is a diagnostic and no strict I claim is made from it.

Use exactly the original positive binary64 Y inputs for a flagged row, save their floathex values, and renormalise them in the fixed high-precision fallback. Do not redraw, reseed, drop a row, clip an observable, or replace a signed gap by its positive part. A nonpositive/nonfinite original Y is a batch error, not a meaningful shape to repair by resampling.

The global maximum of4096 same-input fallback rows is a hard cap. The4097th required row stops the batch with checkpoints and an unresolved/error status; it is not omitted from the sample. High-precision fallback does not turn the finite expectation into a deterministic integral certificate.

## 6. Bounded confidence statements remain diagnostic

For ideal independent P rows, each mathematical level/nonnegative gap X lies in[0,B], B=D/64, so Hoeffding gives

    Pr(|mean X−E X|>=e)<=2exp(−2N e²/B²).

The signed gap has range2B. Bonferroni over a fixed declared list of six levels,eight nonnegative gaps and the signed gap gives simultaneous sampling-model intervals. Use a certified upper B if evaluating this expression numerically. Computed floating means also have arithmetic error; without a certified total arithmetic allowance, label these as nominal bounded-sampling diagnostics rather than unconditional strict expectations. Standard errors or block fluctuations are descriptive, not missing-tail certificates.

A positive sample mean matrix gap, even with a small confidence interval, is not a strict I improvement. A materiality decision may prioritise a route using these diagnostics, but the new upper bound itself still needs a deterministic analytical moment/integration certificate. This theory task executes none of the diagnostic rows.
