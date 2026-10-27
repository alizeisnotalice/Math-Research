# R003: exact Rao–Blackwell truncated Green kernel and certification interfaces

2026-10-27. Independent clone at b2139fc1fe70fd55e7d409031256195c7d10e768. Fixed original64 unnormalised body directions, weight1, common Cauchy time, Lebesgue centred7-cubes and all positive common radii. Only this round's owned directory is written. R001/R002 evidence and other projects are unchanged.

**Conclusion.** The proposed Beta expression is correct, including its epsilon² argument and Gamma(65/2) prefactor. For epsilon>0, H_epsilon is finite and smooth on all R7, including S3. Individual Rao–Blackwell shape terms have uniform spatial amplitude, gradient and Hessian bounds. These permit deterministic interval quadrature plus spatial covering; a raw finite-shape average still does not automatically approximate its expectation with a certified error. A certified lower bound for H_epsilon on the positive face would supply a valid g for the **original G**, and can be combined with R002's original full-G integral budget and the original finite-source/general-L1/all-scale proof. This batch proves no quantitative g and no C7 bound and draws no MC shape.

## 1. Exact model and radial conditioning

Let V have the64 columns v=(1,s2,...,s7), sj=+/-1, so VV^T=64I. Let independent Y_i~Gamma(1/2,1), and S_i=1/(4Y_i). Their Laplace transform is E exp(-lambda S_i)=exp(-sqrt(lambda)), as follows from the elementary integral used in R001/R002. Therefore the Gaussian mixture with covariance2t²Sigma,

    Sigma=sum_i S_i v_i v_i^T,

has characteristic function exp(-t sum_i|v_i·xi|), exactly the original common-time Cauchy process. Its conditional density is

    p_t(x|Sigma)=(4pi t²)^(-7/2)det(Sigma)^(-1/2)
                        exp[-x^TSigma^(-1)x/(4t²)].      (1)

All conditional integrands below are positive. Positive Tonelli, rather than a formal interchange of signed integrals, applies throughout.

Write R=sum_iY_i and P_i=Y_i/R. A direct joint-density change with Jacobian R^63 gives

    R~Gamma(32,1), P~Dirichlet(1/2,...,1/2), independent.

For P in the open simplex define

    A(P)=sum_i v_i v_i^T/P_i, Q_A(x)=x^TA(P)^(-1)x.

Then Sigma=A/(4R). Substitution in (1) gives

    p_t(x|R,P)=pi^(-7/2)det(A)^(-1/2)t^(-7)
                         R^(7/2)exp[-R Q_A(x)/t²].

Performing the Gamma integral in R exactly yields

    p_t(x)=K0 E_P[det(A)^(-1/2)t^(-7)
                       (1+Q_A(x)/t²)^(-71/2)],
    K0=Gamma(71/2)/(pi^(7/2)Gamma(32)).                  (2)

This is Rao–Blackwell conditioning on the same process, not a new shape distribution. In particular the exponent is32+7/2=71/2.

## 2. Truncated-time Beta formula and all constants

Let epsilon>0 and H_epsilon(x)=integral_epsilon^infinity p_t(x)dt. Set

    b=65/2, c=b+3=71/2, d(P)=det(A(P))^(-1/2),
    z=Q/(Q+epsilon²), Q=Q_A(x).

For Q>0 the substitution u=Q/(Q+t²) in the time integral of (2) gives

    integral_epsilon^infinity t^(-7)(1+Q/t²)^(-c)dt
       =(1/(2Q³)) integral_0^z u²(1-u)^(b-1)du.

Indeed t^(-7)dt contributes -u²(1-u)^(-4)du/(2Q³), while the profile supplies (1-u)^(b+3). The exponent remaining is b-1. Since

    B(3,b)=2/[b(b+1)(b+2)],
    (K0/2)B(3,b)=Gamma(b)/(pi^(7/2)Gamma(32))=:KD,

we obtain the exact identity

    H_epsilon(x)=E_P Phi_epsilon(P,x),
    Phi_epsilon(P,x)=KD d(P) Q_A(x)^(-3)
                          I_{Q_A(x)/(Q_A(x)+epsilon²)}(3,65/2),
    KD=Gamma(65/2)/(pi^(7/2)Gamma(32)).                  (3)

Here I_z is the regularized lower incomplete Beta function. With C=Gamma(3)/(4pi^(7/2))=1/(2pi^(7/2)) and D=2Gamma(65/2)/Gamma(32), one has KD=CD. Thus the suggested C*2Gamma(32.5)/Gamma(32) prefactor is correct. There is no missing4 in the argument: Sigma=A/(4R) makes Q_Sigma=4RQ_A, cancelling the4 in (1).

An independent probability derivation checks the same parameter shift. Integrating (1) in t first gives C det(Sigma)^(-1/2)Q_Sigma^(-3) P_gamma(3,Q_Sigma/(4epsilon²)). After substitution it is2C d Q^(-3) E_R[sqrt(R) P_gamma(3,RQ/epsilon²)]. Weighting Gamma32 by sqrt(R) yields Gamma(65/2) times Gamma(65/2)/Gamma32. For an independent U~Gamma3 and T~Gamma(65/2), U/(U+T)~Beta(3,65/2), and U<=TQ/epsilon² is precisely U/(U+T)<=Q/(Q+epsilon²). This reproduces (3).

For integer first parameter3 there is a convenient exact closed form

    I_z(3,b)=1-(1-z)^b[1+bz+b(b+1)z²/2].              (4)

Differentiating the right side gives z²(1-z)^(b-1)/B(3,b), with value0 at z=0 and1 at z=1. It is an identity, not a numerical fit. Formula (4) can lose precision for small z because of cancellation; direct incomplete-Beta interval evaluation or the positive integral in (2) is preferable there.

At Q=0, the apparent Q^(-3) singularity in (3) is removable:

    lim_(Q down0) Q^(-3) I_{Q/(Q+epsilon²)}(3,b)
        =b(b+1)(b+2)/(6epsilon^6).                     (5)

Equivalently

    Phi_epsilon(P,0)=K0 d(P)/(6epsilon^6).

The Gamma recurrence Gamma(b+3)=b(b+1)(b+2)Gamma(b) makes the two expressions identical. Hence (3), with this continuous extension, is valid for every x including0. No shape term is evaluated by dividing0 by0.

## 3. Global continuity and analytic spatial control

The true characteristic exponent is psi(xi)=sum_i|v_i·xi|. From VV^T=64I and l1>=l2,

    psi(xi)>=8|xi|.

The truncated kernel has the absolutely integrable Fourier representation

    H_epsilon(x)=(2pi)^(-7) integral_R7
                  exp(ix·xi-epsilon psi(xi))/psi(xi) dxi.             (6)

Near xi=0 the radial power is r^5; at infinity the exponential dominates all polynomial factors. Therefore H_epsilon belongs to C_infinity(R7), and all spatial derivatives are C0. It is finite at every S3 point where the original G is infinite. It is strictly positive because p_t is strictly positive.

For every unit direction u and integer k>=0, differentiation under (6) gives the deterministic uniform bound

    ||D_u^k H_epsilon||infinity <= D_k(epsilon),
    D_k(epsilon)=Gamma(6+k)/[120pi^4·8^(7+k)epsilon^(6+k)].             (7)

To verify the constant, integrate |xi|^k exp(-8epsilon|xi|)/(8|xi|) in polar coordinates, with |S6|=16pi^3/15 and Fourier factor(2pi)^(-7). In particular

    D0=1/(8^7pi^4epsilon^6),
    D1=6/(8^8pi^4epsilon^7),
    D2=42/(8^9pi^4epsilon^8).

These are upper-bound expressions, not values or lower bounds for H.

The same representation extends holomorphically to the Euclidean tube |Im x|<8epsilon. For |Im x|<=s<8epsilon,

    |H_epsilon(x)|<=D0(epsilon-s/8).

On a face with six complex free coordinates each imaginary part at mosth, it is enough that sqrt(6)h<8epsilon. This supplies a complex-tube bound for a future certified polynomial scheme. A future tensor interpolation implementation must still pay its own dimension-six operator/error constants; none is inherited from n3.

An elementary alternative needs no interpolation theorem: Taylor's formula along a segment gives an Nth-order total-degree Taylor remainder at displacement h bounded by D_(N+1)(epsilon)|h|^(N+1)/(N+1)!. Certified derivatives at the expansion centre plus this bound imply a continuous enclosure. All grid/cover boundaries must be included.

## 4. Uniform bounds on each Rao–Blackwell shape kernel

The cutoff removes the specific obstruction affecting the untruncated full-G shape mean: individual Phi_epsilon(P,·) are uniformly bounded and Lipschitz for all P in the open simplex. This does not itself certify an empirical approximation to their expectation.

First,

    A(P)>=64I, ||A(P)^(-1)||op<=1/64.                  (8)

This follows from P_i<=1 and the Gram identity. A useful stronger **determinant** bound is

    det A(P)>=64^14, d(P)<=64^(-7).                    (9)

For completeness, let alpha_J=det(V_J)^2/64^7 for7-column subsets. Cauchy-Binet gives sum alpha=1 and frame leverage gives sum_(J contains i)alpha_J=7/64 (differentiate det(sum t_i v_iv_i^T) at t_i=1). Weighted AM-GM then gives

    det A>=64^7 product_i P_i^(-7/64)>=64^14,

where the last step uses product_iP_i<=64^(-64). Equality holds at P_i=1/64, where A=4096I. No enumeration of all7-subsets is required. The determinant lower bound (9) does not claim every eigenvalue is at least4096; inverse estimates below use only (8).

Write again c=71/2 and define

    kappa_c=(2c+1)^(-1/2)[(2c+1)/(2c+2)]^(c+1),
    theta_c=[(c+1)/(c+2)]^(c+2).

Uniformly for all P and x,

    0<Phi_epsilon(P,x)<=B_epsilon,
    B_epsilon=K0/[6·64^7 epsilon^6],
    ||grad_x Phi_epsilon(P,x)||<=L_epsilon,
    L_epsilon=K0 c kappa_c/[28·64^7 epsilon^7],
    ||D_x²Phi_epsilon(P,x)||op<=J_epsilon,
    J_epsilon=K0(4c theta_c+2c)/[512·64^7 epsilon^8].    (10)

The amplitude follows directly from (2) and (9). To prove the gradient, for any unit u put r=sqrt(Q)/t. By (8),

    |u^TA^(-1)x|<=sqrt(u^TA^(-1)u)sqrt(Q)<=sqrt(Q)/8.

Differentiating (1+Q/t²)^(-c) therefore bounds its directional derivative by(c/(4t))r(1+r²)^(-c-1). The supremum of the last r-profile is kappa_c, attained at r²=1/(2c+1). Integrate t^(-8) from epsilon to infinity and use (9).

For the Hessian in unit direction u, the two differentiated terms have absolute sum at most

    1/(64t²) [4c(c+1)r²(1+r²)^(-c-2)
                         +2c(1+r²)^(-c-1)].

The first r-profile has supremum theta_c/(c+1), attained at r²=1/(c+1), and the second is at most1. Integrating t^(-9) gives the stated J_epsilon. Directional quadratic-form control gives the Hessian operator norm because the Hessian is real symmetric.

Bounds (10) justify differentiation through the P expectation and give the same B,L,J bounds for H_epsilon. For a spatial covering one can use the smaller of L_epsilon and D1(epsilon), and similarly the Hessian bounds. The finite check records them only at the registered epsilon=1/10; approximately B=0.00014594534, L=0.00079085577, J=0.02095704149. These are deliberately uniform bounds and may be too coarse for a useful six-dimensional cover.

For interval implementation it is helpful to retain the positive t-integral. The scalar profile F_epsilon(Q)=Q^(-3)I_(Q/(Q+epsilon²))(3,b), extended at0, is decreasing, with alternating nonnegative derivatives: (-1)^k F_epsilon^(k)(Q)>0. Differentiate (2), whose kth Q derivative multiplies by (-1)^k(c)_k t^(-2k). At Q=0 the exact corresponding Phi derivative is

    partial_Q^k Phi(P,0)
      =(-1)^k K0 d(P)(c)_k /
                    [(6+2k)epsilon^(6+2k)].           (11)

Thus Q intervals can be evaluated monotonically, including0, instead of subtracting nearly equal terms in (4). The x derivatives are

    grad Phi=2KD d F'(Q) A^(-1)x,
    Hess Phi=2KD d F'(Q) A^(-1)
           +4KD d F''(Q)(A^(-1)x)(A^(-1)x)^T.

An outward interval LDL/Cholesky SPD certificate, determinant and solve, together with regularized incomplete Beta or (11), can enclose these. Floating SPD/conditioning diagnostics alone are not such a certificate; no clipping or removal of ill-conditioned shapes is allowed.

## 5. Deterministic shape-integration error on a compact bulk

Let K_delta={P: P_i>=delta for every i}, with0<delta<1/64. Use simplex chart coordinates p1,...,p63 and p64=1-sum_(j<=63)p_j. This is convex. On this set

    partial_(p_j)A=-v_jv_j^T/p_j²+v_64v_64^T/p_64²,
    ||partial_(p_j)A||op<=14/delta².                    (12)

For a symmetric matrix perturbation E, differentiating the determinant/profile in (2), using (8), gives

    |D_A Phi[E]| <=(39/64)B_epsilon||E||op.             (13)

Indeed the logarithmic determinant term is bounded by7||E||/(2·64). The quadratic-form variation satisfies |D_A Q[E]|=|x^TA^(-1)EA^(-1)x|<=||E||Q/64. Hence the profile term contributes at most c||E||/64, since(Q/t²)/(1+Q/t²)<=1. Their sum is(7/2+c)/64=39/64. Multiplying by the undifferentiated positive upper envelope and integrating proves (13).

Combining (12)-(13), uniformly in x and every j<=63,

    |partial_(p_j)Phi|<=L_P(delta),
    L_P(delta)=(273/32)B_epsilon/delta².                (14)

Thus |Phi(P,x)-Phi(P',x)|<=L_P(delta)sum_(j<=63)|P_j-P'_j| for P,P' in K_delta. This estimate includes the chart-dependent P64 term, which must not be omitted.

Suppose disjoint measurable cells C_j in K_delta have exact masses w_j under the Dirichlet law, anchors P_j and chart l1 radii r_j. Then uniformly in x,

    |integral_Kdelta Phi(P,x)dmu -sum_j w_jPhi(P_j,x)|
       <=L_P(delta)sum_j w_jr_j                       (15)

if the cells cover K_delta up to null sets. If they cover only a subset, (15) is valid for that subset and positivity allows a **lower** certificate without paying the missing positive mass. For a two-sided full expectation enclosure use the tail bound

    mu(K_delta^c)<=min(1,64 I_delta(1/2,63/2))=:tau_delta,
    integral_Kdelta^c Phi dmu<=B_epsilon tau_delta.     (16)

Each P_i has Beta(1/2,63/2) marginal, so (16) is just a union bound; it does not multiply dependent-coordinate probabilities. If weights themselves are approximated, add B_epsilon times the l1 weight error. For very small delta the factor delta^(-2) can make (14) impractical. The existence of a deterministic formula does not imply that a63-dimensional partition will have a feasible cell count.

An explicit source of exactly certifiable cell masses is independent stick-breaking:

    U_j=P_j/(1-sum_(i<j)P_i), j=1,...,63,
    U_j independently~Beta(1/2,(64-j)/2),
    P_j=U_j product_(i<j)(1-U_i),
    P64=product_(i<=63)(1-U_i).

This follows by the same Gamma/Dirichlet density change. Rational boxes in U have masses equal to products of63 Beta-CDF differences, which Arb can enclose outward. Boxes must be disjoint and checked to map into K_delta; no independence of the P_i is assumed. Their images need not be convex: their anchors and points lie in the convex K_delta, so (14) still applies along the connecting segment.

For a face patch X_nu of Euclidean radius rho_nu about x_nu, point enclosures and lower mass bounds produce the fully explicit certificate

    H_epsilon(x)>=sum_j w_j^-
       max(0, Phi(P_j,x_nu)^- -L_P(delta)r_j
                                      -L_epsilon rho_nu), x in X_nu. (17)

Here cells need only be disjoint selected positive-mass cells, not a full shape cover. The positive part ensures lower mass endpoints are multiplied by nonnegative contributions. A certified Hessian/Taylor enclosure can replace the Lipschitz losses. Formula (17) pays shape-cell error, arithmetic error and spatial coverage separately. A sample mean with equal weights supplies none of its mass-cell statements automatically.

## 6. From a face H lower certificate to the original g/I transfer

For fixed P, Phi_epsilon(P,x) is decreasing in Q_A(x). Consequently it is radially nonincreasing in x, so H_epsilon(r y)>=H_epsilon(y) for0<=r<=1. The true H_epsilon also has the full cube sign/permutation symmetry because the original common-time process does; an arbitrary finite selected-shape lower surrogate need not have that symmetry.

Therefore a continuous certificate

    H_epsilon(1,y)>=g0>0 for every y in [0,1]^6

extends by the true symmetry to every boundary face and then by radial monotonicity to the entire cube. Since H_epsilon<=G pointwise by positive time truncation, it implies G>=g0 on the cube. Alternatively one may use the homogeneity of **G**, G(ry)=r^(-6)G(y); one must not incorrectly call H_epsilon homogeneous. Its exact scaling is

    H_epsilon(r x)=r^(-6)H_(epsilon/r)(x), r>0.

A finite list of evaluated face points is insufficient. With a certified point expectation h_nu^- and a covering of radius rho_nu, H>=h_nu^- -min(L_epsilon,D1)rho_nu suffices on the patch. But the h_nu values must include the deterministic P-expectation error, not only the interval error of a finite-shape mean.

Once such g0 exists, keep the original G in the majorant construction:

    c=g0/7,
    K=(G-c)_+/[128(g0-c)].

R001 proves its whole-space finite source, positive probability measure and general-L1/all-r comparison for the original generator. R002 proves the full original-G integral bound

    I=integral_[0,1]^6 G(1,y)^(7/6)dy
       <R_I=178444228577/20000000000000000.

They combine legitimately with the H-derived g0 to give the conditional bound

    C7<=7^(7/6)R_I/g0^(7/6).

No short-time norm must be added to R_I: R_I already covers the complete original G and all positive times. If a future method instead certifies an integral of H only, it must separately pay the R001 whole-face short-time L^(7/6) error.

It is **not** valid to replace G by H_epsilon in the positive-part source argument without a new proof. In fact

    L H_epsilon=-p_epsilon

on the whole space, by integrating the semigroup derivative from epsilon to infinity (or directly from (6)). Its negative source is a spread-out density, not delta0. Its positive-part source is not the original mu-delta0 identity. The legal role here is a lower certificate for original G, while the original G supplies its own finite source, mass identity, common-radius scaling, general L1 limit and all-scale comparison.

Since H_epsilon is strictly positive continuous on the compact cube, its actual minimum is positive. This qualitative existence statement is not a numerically useful g certificate and is not entered into the bounds ledger. No value of that minimum has been certified in this batch.

## 7. A concrete next finite computation, and what it would accept

The next useful computational object is an **outward interval evaluator of Phi_epsilon and its derivatives**, not another optimization of the full-G cache. A separately registered deterministic pilot can use epsilon=1/10, rational shapes and face patches, with the following finite acceptance problem:

1. Implement (3), its Q=0 branch and the positive derivative formulas. Use a fixed rational shape set containing uniform P and one shape biased toward each of the8 Hadamard blocks (for a chosen block, P_i=1/16 on its8 columns and1/112 on the other56). These9 rational shapes require no random draw. Evaluate at a fixed rational face point and at0, retaining interval SPD/solve certificates and comparing Beta against the positive time integral or derivative identity. The large shape amplitude is finite even at S3.
2. Construct **one** rational stick-breaking box around a chosen shape, with its63 Beta-CDF mass factors, a checked delta and chart radius. For **one** closed face patch, certify (17) or a Hessian improvement. Return its lower contribution, all weight/kernel intervals and the cell/patch geometry. This validates the complete deterministic lower-certificate pipeline. It need not be a useful global g and must not be presented as one.
3. Before any bulk expansion, compare the required global target to the provable pilot losses. For the generic128 threshold the conditional required floor is g_req=7(R_I/128)^(6/7). Only if the certified shape-mass and spatial losses leave enough margin is a new finite full-face cover justified; otherwise record that this specific cell construction is too weak and stop. Very small box mass or a huge delta^(-2) loss is a useful explicit rejection, not permission to grow an unbounded63-dimensional partition.

This is a reproducible finite interval task with pass/fail conditions. It produces no expectation error from MC standard errors and cannot be replaced by interval arithmetic on the existing8192 sample average. Existing frozen shapes could still be used for diagnostic H_epsilon values after separate counting, but that alone would not satisfy the deterministic weight/error condition. No such pilot, old-cache evaluation or full-face run is executed in the present theory batch.

## 8. Actual checks, provenance and scope

`registration.md` was written before the finite check. The first actual `mincheck.py` run exited0 in0.12568545900285244 script seconds, with no failed assertion or repair. It exactly reconstructed Gram64I and uniform-P A=4096I, checked the closed Beta derivative, endpoint values, Q=0 coefficient and Gamma recurrences. At512-bit Arb it compared (4) with direct regularized Beta at exactlyQ/epsilon²=0,1/64,1,64; the interval differences contain0. These numerical overlaps diagnose the formula implementation; the identity itself is proved analytically in Sections1-2.

At the registered epsilon=1/10 it enclosed the B,L,J and Fourier-order0/1/2 upper-bound expressions. These are deterministic **bound expressions**, not sampled H values, expectation intervals, g lower values or C7 results. All source hashes were unchanged before/after. No MC sample, full-G evaluation, MATLAB call, commit or push occurred.

Reproduce from this independent clone:

    /Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python -B cube-dim-7-bounds/upper/round003/truncated_green/mincheck.py

`results.json` includes exact outward rational endpoints and versions; `sources.json` locks the read-only sources; stdout/stderr and `receipt.json` record execution; `manifest.json` lists the final files. The theoretical derivation and uniform error interfaces are complete for this registered batch. The missing expectation quadrature, full-face covering and useful g remain future work. Stop for main-agent review.
