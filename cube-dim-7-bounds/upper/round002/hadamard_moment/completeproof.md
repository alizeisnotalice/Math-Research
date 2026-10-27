# Round002: Hadamard-block determinant moment and complete I budget

2026-10-27. Fixed n=7, equal-weight64 unnormalised body directions, the original common-time Cauchy process, Lebesgue measure, centred axis-parallel cubes. This independent new clone writes only this directory. Round001 is read-only. No Monte Carlo, old Monte Carlo baseline, Green node, MATLAB, agent expansion, model change, commit or push occurs here.

**Result.** The64 directions split into8 disjoint tight7x8 Hadamard blocks. A self-contained SPD determinant inequality and independent Gamma moments give the rigorous analytic endpoint

    M=E det(Sigma)^(-1/14) <= M_H,
    M_H=8^(-1/2) [2^(-1/28) Gamma(65/16)/Gamma(57/14)
                      ·(Gamma(57/112)/Gamma(1/2))^8]^8.

The512-bit enclosure proves M_H<M_Gamma, with M_H approximately0.08730273429780108714944687916 and M_Gamma approximately0.09456236373702571242462782592. Consequently the complete time-and-boundary integral satisfies

    I <= I_H := M_H^(7/6)/(105·2^(25/6)·pi^(13/12))
      < 178444228577/20000000000000000
       =0.00000892221142885.

The last rational endpoint is also strictly below the former analytic I_Gamma expression. This certifies an improved **complete I upper budget**, not an estimate of actual M or actual I from below. No globally useful g lower certificate is provided; **there is no new C7 upper bound**. The registered successful batch stops here without further exponent variants.

## 1. Original model and Gaussian mixture

Write v=(1,s2,...,s7), sj in {-1,1}, and V for these64 columns. The original process has characteristic function

    exp(-t sum_v |v·xi|).

Let Y_i be iid Gamma(1/2,1), S_i=1/(4Y_i), and Sigma=sum_i S_i v_i v_i^T. All S_i are strictly positive almost surely. The Laplace transform of S is exp(-sqrt(lambda)). For completeness, put

    J(c)=integral_0^infinity exp(-u²-c²/u²)du, c>0.

Differentiation and the change w=c/u give J'(c)=-2J(c); dominated convergence at c=0 gives J(0)=sqrt(pi)/2. Therefore J(c)=(sqrt(pi)/2)exp(-2c), and the substitution Y=u² yields

    E exp(-lambda/(4Y))=exp(-sqrt(lambda)).

All differentiations are justified on compact positive c intervals by the exponential decay at0 and infinity. This establishes the mixture without changing the time convention.

Given Sigma, use the centred Gaussian of covariance2t²Sigma. Its characteristic function is exp(-t²xi^TSigma xi). Averaging over the independent S_i gives exactly the original common-time characteristic function. Its conditional density is

    (4pi t²)^(-7/2)det(Sigma)^(-1/2)
       exp[-x^TSigma^(-1)x/(4t²)].

It is an actual positive Gaussian density; no signed fan or separate shape process is introduced. Define G(x)=integral_0^infinity p_t(x)dt as an extended nonnegative function. With Q_Sigma(x)=x^TSigma^(-1)x, positive Tonelli and u=Q/(4t²) give, for x!=0,

    G(x)=C E[det(Sigma)^(-1/2) Q_Sigma(x)^(-3)],
    C=Gamma(3)/(4pi^(7/2))=1/(2pi^(7/2)).                 (1)

Any infinite values on singular direction spans are allowed in this identity. The whole-boundary argument below proves the relevant integral finite without discarding their neighbourhoods.

## 2. Exact partition into8 Hadamard blocks

Index the columns of a Sylvester8x8 Hadamard matrix by u in F2^3 and its rows by r in F2^3:

    H_(r,u)=(-1)^(r·u).

Retain rows r=0,1,...,6, deleting row7. The first retained row is all1, and the remaining rows are six distinct nonzero characters. Its8 direction columns can be identified with the binary linear code

    C={ (r·u)_(r=1,...,6): u in F2^3 } subset F2^6.

This map is injective since rows1,2,4 contain the three coordinate characters. Thus C is a3-dimensional linear code with8 words; F2^6/C has8 cosets. The binary masks used by the checker are just an exact representation of these vectors, not a search parameter.

Each coset a+C corresponds to multiplying the six free-coordinate rows by the signs (-1)^(a_r), leaving the first row unchanged. Hence its7x8 matrix H_a is a row-sign transform and a column permutation of the retained Hadamard matrix, and

    H_a H_a^T=8 I7.                                      (2)

The8 cosets are disjoint and together contain every64 body direction exactly once. Each H_a spans R7.

Every7-column minor of H_a has absolute determinant512. An analytical check follows from H8H8^T=8I8: |det H8|=8^4=4096 and H8^(-1)=H8^T/8, so each cofactor has absolute4096/8=512. Deleting row7 and any column gives precisely the canonical7-minors. Row signs and column permutations preserve the absolute value. The checker independently verifies all8 canonical minors and56 corresponding minors in the other7 cosets,64 total, as exact integers; it also checks (2) for all8 blocks and the exact64-column partition.

Let Sigma_g=sum_{i in block g} S_i v_i v_i^T. Each is SPD almost surely by positive weights and (2), they are independent, and Sigma=sum_(g=1)^8 Sigma_g. A row-sign transform preserves determinant and iid S_i are insensitive to a column permutation, so their determinant laws are identical.

## 3. Self-contained determinant Minkowski inequality for SPD matrices

For SPD n x n matrices A,B, let lambda_i>0 be the eigenvalues of A^(-1/2)BA^(-1/2) and rho=(product_i lambda_i)^(1/n). The kth elementary symmetric sum obeys

    e_k(lambda)>=binom(n,k)rho^k, 0<=k<=n.               (3)

For1<=k<=n-1, apply ordinary AM-GM to its binom(n,k) positive summands. Each lambda_i occurs in exactly binom(n-1,k-1) such products, so their geometric mean is (product lambda)^(k/n)=rho^k. The endpoint cases are identities. Summing (3) gives

    det(I+A^(-1/2)BA^(-1/2))
       =product_i(1+lambda_i)>=(1+rho)^n.

Multiplying by det A and taking positive nth roots proves

    det(A+B)^(1/n)>=det A^(1/n)+det B^(1/n).             (4)

This proof supplies both the inequality and the needed hypotheses, so no unverified literature quotation is necessary. Iterating (4) on the8 SPD blocks yields

    det(Sigma)^(1/7)>=sum_g det(Sigma_g)^(1/7).          (5)

All matrices here are SPD almost surely; a semidefinite limiting variant is not needed.

Put W_g=det(Sigma_g)^(1/7)>0. Since x->x^(-1/2) is decreasing, (5) gives

    det(Sigma)^(-1/14)<=(sum_g W_g)^(-1/2).

Ordinary AM-GM has sum_g W_g>=8(product_g W_g)^(1/8). Inverting and taking the square root therefore gives, with the direction explicit,

    det(Sigma)^(-1/14)
       <=8^(-1/2) product_g W_g^(-1/16)
        =8^(-1/2) product_g det(Sigma_g)^(-1/112).        (6)

Thus independence of the disjoint groups implies

    M<=M_H=8^(-1/2)[E det(Sigma_g)^(-1/112)]^8.          (7)

No independence of matrices derived from overlapping groups is assumed. The group covariance matrices are not claimed to be scalar multiples of I for random weights; only the deterministic unweighted frame is tight.

## 4. Exact negative moment of one block

Fix a block of8 variables. Cauchy-Binet and all8 squared minors equal512²=2^18 give

    det(Sigma_g)=2^18 sum_(omit j) product_(i!=j) S_i.

Substituting S_i=1/(4Y_i) yields

    det(Sigma_g)=16 (sum_i Y_i)/(product_i Y_i),          (8)

because2^18/4^7=16. Every term and Jacobian factor is positive. This identity was checked exactly at the coefficient/degree level; it is not fitted to random matrices.

Let T=sum_iY_i and P_i=Y_i/T. The Gamma joint density under the positive transformation y_i=T p_i, whose Jacobian is T^7, factors into

    T~Gamma(4,1), P~Dirichlet(1/2,...,1/2), independent.

The density factorisation directly proves this fact. For every a>0, (8) gives

    det(Sigma_g)^(-a)=2^(-4a) T^(7a) product_i P_i^a.

The exponent7a is8a-a, so the sum in (8) has been paid exactly once. The Gamma integral and the Dirichlet beta integral are

    E T^(7a)=Gamma(4+7a)/Gamma(4),
    E product_iP_i^a=Gamma(4)/Gamma(4+8a)
                        ·[Gamma(1/2+a)/Gamma(1/2)]^8.

For clarity the beta integral itself follows by applying the transformation y_i=T p_i to the product of the8 gamma integrals with shapes1/2+a, and then dividing by the normalising Dirichlet integral. All parameters are positive and both integrals finite. Therefore

    E det(Sigma_g)^(-a)
      =2^(-4a) Gamma(4+7a)/Gamma(4+8a)
                       ·[Gamma(1/2+a)/Gamma(1/2)]^8.    (9)

Inserting a=1/112 gives4+7a=65/16,4+8a=57/14 and1/2+a=57/112. This proves the endpoint stated at the start. Equivalently,

    M_H=2^(-25/14)[Gamma(65/16)/Gamma(57/14)]^8
                        ·[Gamma(57/112)/Gamma(1/2)]^64.

The negative moment formula holds for every a>0, not just the registered exponent. The prescribed numerical batch evaluates only a=1/112.

## 5. Independent analytical reconstruction of the comparison endpoint

For comparison, reconstruct the earlier analytical M_Gamma expression. Let m=64, c_J=det(V_J)^2 for7-column subsets J. Cauchy-Binet at the unweighted tight frame VV^T=mI gives sum_J c_J=m^7. Define alpha_J=c_J/m^7, nonnegative with sum1. Differentiation of det(sum_i t_i v_i v_i^T) at all t_i=1 gives

    sum_(J contains i)alpha_J
       =tr((mI)^(-1)v_i v_i^T)=7/m.

Weighted AM-GM now implies

    det(Sigma)=m^7 sum_J alpha_J product_(i in J)S_i
       >=m^7 product_i S_i^(7/m).

Raising to -1/14 and using independence gives

    M<=m^(-1/2)[E S^(-1/(2m))]^m
      =(1/4)[Gamma(65/128)/Gamma(1/2)]^64
      =:M_Gamma.                                       (10)

Indeed E S^(-1/128)=4^(1/128)Gamma(65/128)/Gamma(1/2). The outer factor is m^(-1/2)4^(64/128)=1/4. This is an analytical endpoint, not an empirical baseline. Equations (7) and (10) are both upper bounds on the same unknown actual M; comparing them cannot give a lower bound on M or on their gaps to M.

The actual512-bit interval calculation proves M_Gamma-M_H>0. Displayed approximations are

    M_H=0.08730273429780108714944687916313938...,
    M_Gamma=0.09456236373702571242462782591738934...,
    M_H/M_Gamma=0.92322918810052829431138774002019656....

Exact outward rational interval endpoints are saved in `results.json`; decimal displays are not the strict comparison certificate.

## 6. Complete boundary and time integral, with no facewise shape symmetry

Let K=[-1,1]^7, q=7/6, and dσ be ordinary6-dimensional face measure on its full boundary. Fix any SPD Sigma, not necessarily symmetric under coordinate permutations. Cubical rays z=rx, x on the boundary, have Jacobian r^6 dr dσ. Hence

    integral_boundary Q_Sigma(x)^(-7/2)dσ
      =[2/Gamma(7/2)] integral_R7 exp(-Q_Sigma(z))dz
      =H sqrt(det Sigma),
    H=2pi^(7/2)/Gamma(7/2)=16pi^3/15.                   (11)

The identity follows by integrating exp(-r²Q(x))r^6 in r; the Gaussian integral on R7 is pi^(7/2)sqrt(det Sigma). Boundaries between face cones have zero measure. No single positive face has been replaced by its symmetry average at a fixed nonsymmetric shape.

For the shape kernel in (1), since3q=7/2, (11) gives

    ||det(Sigma)^(-1/2) Q_Sigma^(-3)||_(Lq(boundary K))
       =H^(6/7)det(Sigma)^(-1/14).

Minkowski's integral inequality therefore yields

    ||G||_(Lq(boundary K))<=C H^(6/7) M.                (12)

There is no circular finiteness assumption. First restrict the expectation to increasing events

    E_N={N^(-1)I<=Sigma<=NI}.

The restricted kernels are bounded on the compact boundary, so Minkowski applies there. Since the moments of det(Sigma)^(-1/14) are bounded by the finite (7), monotone convergence as N increases proves (12) and G in Lq(boundary K). Positive Tonelli in (1) already includes every time0<t<infinity. Thus all angular singularities and their neighbourhoods are included, not merely a finite cutoff or an a.e. deletion of their neighbourhoods.

Only **after taking expectation**, use the original frame's full coordinate sign/permutation symmetry. Signed permutations map its64 unoriented body lines to themselves and preserve the iid weight law. Consequently, with

    I=integral_[0,1]^6 G(1,y)^q dy,

all14 faces and their64 orthants have the same integral:

    integral_boundary G^q dσ=14·64 I=896 I.

Raise (12) to q and use C=1/(2pi^(7/2)), H=16pi^3/15:

    I<=C^q H M^q/896
      =M^(7/6)/(105·2^(25/6)·pi^(13/12)).              (13)

The power of pi is -49/12+3=-13/12; the factor is16/(896·15)·2^(-7/6)=1/(105·2^(25/6)). Formula (13) has been rederived here rather than inferred from the Round001 status label.

Combining (7) and (13) gives the complete analytic I_H expression. The checker evaluates (13) both at M_H and at the comparison expression M_Gamma. It proves

    I_H < R_I :=178444228577/20000000000000000 < I_Gamma,
    I_H approximately8.92221142884962178793463878578e-6,
    I_Gamma approximately9.79365339029564320881646517954e-6.

The rational endpoint is an explicit outward upper bound obtained from the actual512-bit interval, rounded upward to a fixed10^(-18) grid. The saved exact dyadic upper endpoint is finer. In particular the actual mathematical I obeys I<=I_H<R_I. None of the lower endpoints in `results.json` is asserted to lower-bound the unknown true I or M.

## 7. Meaning for the upper-bound programme and stopping boundary

For a separately proved global g0<=essinf_K G, the accepted original-process finite-source interface gives

    C7<=7^(7/6) I0/g0^(7/6).

This work supplies I0=R_I with a strict underlying expression bound, including all positive times and all faces. It supplies no useful certified g0. A finite point, an empirical minimum or an unproved symmetry reduction would not close that missing obligation. Thus no new numerical C7 upper bound is claimed, and the current C7 threshold is not updated by this directory.

The successful prescribed AM-GM expression removes any reason to open additional Jensen/Holder weight variants in this batch. There were no Monte Carlo shapes, new Green samples, baseline science runs, quadrature, model optimization or adaptive extension. The determinant inequalities, exact moments and full-boundary identity are the continuum proof; the finite computation only certifies the partition, integer constants and strict scalar interval comparisons.

## 8. Actual finite run and reproducibility

The finite parameters and stop conditions were written before execution in `registration.md`. The first and only registered check run completed exit0 in0.058543666993500665 script seconds, with no mathematical assertion failure or implementation repair. Its64 minors mean8 canonical plus56 additional minors, matching the intended8-coset verification rather than accidentally omitting the canonical block. The source hashes were taken before and after and remained identical.

Run from this independent clone:

    /Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python -B cube-dim-7-bounds/upper/round002/hadamard_moment/mincheck.py

`results.json` stores all8 integer blocks,64 signed minors, binary cosets, every scalar interval with exact rational endpoints, the explicit I upper rational, Python/SymPy/python-flint versions,512-bit precision, code/source hashes and measured runtime. `run.stdout` is its short actual output and `run.stderr` is empty. `sources.json` locks only the read-only local analytical sources; no old source or large historical evidence is copied here. `receipt.json` and `manifest.json` record the run and final whitelist.

The proof of the SPD determinant theorem is self-contained in Section3, so its validity does not depend on a secondary theorem citation or on a missing external page. SymPy/FLINT/Arb are shared libraries; this is not a second arithmetic kernel. All new files remain local for main-agent integration and review. This round stops pending that review.
