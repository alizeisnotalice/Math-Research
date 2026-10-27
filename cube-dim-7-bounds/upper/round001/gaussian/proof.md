# Seven-dimensional full-boundary Gaussian-mixture moment interface

This is a new-dimension proof and fresh finite probe of the existing n=6 R123 method, not a claim of a new mechanism. All analytic statements below concern the original equal-weight body64, common-time Cauchy model. The only inherited method source read was `cube-dim-6-bounds/upper/round123/shape_Minkowski_moment/theorem.md`; its historical numerical log is not new evidence here. No fan, altered generator, facewise single-shape symmetry, or sampled global minimum is used.

## 1. Exact model and Gaussian representation

Let n=7, m=64, and let V have columns v=(1,s_2,...,s_7), s_i in {-1,1}. Thus VV^T=m Id. The Cauchy semigroup has characteristic function exp(-t sum_v |v dot xi|), corresponding to independent standard Cauchy variables with the same time t and unit weights. Its potential is G(x)=integral_0^infinity p_t(x) dt (extended positive values allowed).

Take iid S_v>0 with E exp(-lambda S)=exp(-sqrt(lambda)). Given Sigma=sum_v S_v vv^T, the centered Gaussian with covariance 2t^2 Sigma has characteristic function exp(-t^2 xi^T Sigma xi). Averaging gives exactly the displayed original semigroup. Its density is (4pi t^2)^(-7/2) det(Sigma)^(-1/2) exp(-Q/(4t^2)), where Q=x^T Sigma^(-1)x. Tonelli and u=Q/(4t^2) give

G(x)=C E[det(Sigma)^(-1/2) Q^(-3)],  C=Gamma(3)/(4pi^(7/2))=1/(2pi^(7/2)).     (1)

For general n>1 the same calculation gives Gamma((n-1)/2)/(4pi^(n/2)), so no new time convention has entered.

## 2. Whole-boundary identity and Minkowski

Use ordinary six-dimensional face measure on boundary K, K=[-1,1]^7; intersections of faces are null. For any SPD Sigma, cubical ray coordinates z=rx, x in boundary K, have Jacobian r^6 dr dσ. Hence

integral_boundary Q(x)^(-7/2) dσ = [2/Gamma(7/2)] integral_R7 exp(-Q(z)) dz
 = H sqrt(det Sigma),  H=2pi^(7/2)/Gamma(7/2)=16pi^3/15.                 (2)

This holds for every shape, without coordinate symmetry of that shape. Put q=7/6 and M=E det(Sigma)^(-1/14). Lq Minkowski on the **entire** boundary gives

||G||_q <= C H^(1/q) M.                                               (3)

Indeed the shape Lq norm is H^(1/q) det(Sigma)^(1/(2q)); multiplication by det(Sigma)^(-1/2) leaves exponent -1/14. To avoid assuming finiteness beforehand, restrict the expectation to increasing events {N^(-1)Id <= Sigma <= N Id}, apply Minkowski there, and then monotone convergence to the q-th powers. Section 3 proves M finite.

The original G, after expectation, has coordinate sign and permutation symmetry because these permute the 64 unoriented body lines with equal weights. Define I=integral_[0,1]^6 G(1,y)^q dy. Only at this stage use integral_boundary G^q=7*2^7 I=896I. Thus

I <= I_M(M):= C^q H M^q/896
 = M^(7/6)/(105*2^(25/6)*pi^(13/12)).                                  (4)

Equation (4) includes all positive times, the complete boundary, and its angular singularities. An empirical M is not a rigorous expectation upper endpoint.

## 3. Exact Gamma/Dirichlet radial elimination and analytic upper endpoint

Y_v=1/(4S_v) are iid Gamma(1/2,1). With R=sum Y_v and P_v=Y_v/R, R~Gamma(32,1) is independent of P~Dirichlet(1/2,...,1/2). Set A(P)=sum_v vv^T/P_v and Z(P)=det(A)^(-1/14). Then Sigma=A/(4R), and

M=D E_P Z,  D=2 Gamma(65/2)/Gamma(32).                                 (5)

For point evaluation put T_x(P)=det(A)^(-1/2)(x^T A^(-1)x)^(-3). The exact same radial elimination gives

G(x)=C D E_P T_x(P),                                                   (6)

as an extended positive expectation. No independent R is sampled. Equations (5) and (6) are Rao--Blackwellization, not a change of model.

Here is a direct deterministic envelope. For seven-column subsets J, write c_J=det(V_J)^2. Cauchy--Binet gives sum_J c_J=m^7 and alpha_J=c_J/m^7 are nonnegative weights with sum 1. Differentiating det(sum_i t_i v_i v_i^T) at all t_i=1 shows sum_{J contains i}alpha_J=7/m (matrix determinant derivative, since |v_i|^2=7 and VV^T=mId). Weighted AM--GM therefore gives

det A=m^7 sum_J alpha_J product_{i in J} P_i^(-1)
 >=m^7 product_i P_i^(-7/m).

Consequently

0<Z<=m^(-1/2) product_i P_i^(1/(2m))<=1/m=1/64.                         (7)

Equality holds at uniform P, where A=m^2 Id. The Dirichlet beta integral yields

E product_i P_i^(1/128)
 = Gamma(32)/Gamma(65/2) * [Gamma(65/128)/Gamma(1/2)]^64.

Combining with (5),

M <= M_Gamma := (1/4)[Gamma(65/128)/Gamma(1/2)]^64 < infinity.            (8)

This explicit analytic endpoint inserted into (4) is a rigorous complete I upper bound. The certificate checker evaluates its enclosing Arb interval from exact rational inputs. It does not prove a global g floor. The simpler endpoint M<=D/64 is also valid. As an independent lower diagnostic, Hadamard and all diagonal entries Sigma_ii=sum S_v give det Sigma<=(sum S_v)^7. The sum has law m^2 S, so M>=E S^(-1/2)/m=1/(32sqrt(pi)). This is not an operator lower bound.

## 4. Positive-part mass and the general L1 interface

G is homogeneous of degree -6. For c>0 and g=G(x) on the boundary, cubical rays and Tonelli give

integral_0^infinity (r^(-6)g-c)_+ r^6 dr=(6/7)c^(-1/6)g^(7/6).

This is valid in the extended sense and (4) makes its boundary integral finite. Thus

integral_R7 (G-c)_+ = (6/7)c^(-1/6)*896I.                               (9)

If a separate proof supplies g_L>0 with G>=g_L on K, then 0<c<g_L gives the nonnegative majorant kernel (G-c)_+/[128(g_L-c)] >=1_K/128. Its complete mass is 6 c^(-1/6)I/(g_L-c), minimized at c=g_L/7. The resulting interface is

U(g_L,M_U)=7^(7/6) I_M(M_U)/g_L^(7/6).                                 (10)

This only replaces the I estimate in the accepted original-generator positive-source transfer; it does not independently reprove that transfer. With that transfer and a certified g_L, (10) applies to the Lebesgue centered axis-parallel cube maximal operator on general L1 and all positive common radii. No finite-point G estimate supplies such a g_L. Without it, no new C7 bound follows here.

## 5. Exact rank audit and the fixed eight candidate points

The initially requested anchors a_k=(1,1 repeated k times,0 repeated 6-k times), k=0,...,6, are audited without Monte Carlo after the parent's pre-run revision. For k<6, a_k=(v+w)/2 with v,w equal to 1 on its 1-coordinates and opposite on its 0-coordinates. Thus they lie in a two-body span (and k=6 is a single body vector). Moreover G(a_k)=+infinity: remove these one or two Cauchy variables; the remaining 63 or 62 body variables still span R7, so their time-one density is continuous and bounded below near zero. Keeping the selected Cauchy coordinates within O(t) of 1 or 1/2 has probability at least const*t^2 per selected nonzero coefficient. In the convolution density this gives p_t(a_k)>=const*t^(2r-7), r=1 or 2, for small t. Both powers integrate to infinity. To justify the remainder density statement, select seven independent remaining columns; their independent Cauchy variables have an everywhere positive continuous joint density, and convolution with the remaining variables preserves positivity, continuity and a local positive lower bound. No anchor is a finite g diagnostic.

The eight fixed rational face points are x_j=(1,(j+1)/(j+7),...,(j+6)/(j+7)), j=0,...,7. Each has seven distinct positive absolute coordinate values. A vector in the span of three body lines has at most four different absolute coordinate values: the eight signed sums ±a±b±c form at most four pairs under negation. Therefore every x_j avoids the union S3 of all at-most-three-body spans. This exact fact alone is not asserted here to prove all pointwise integration estimates or a global floor. Their same-shape estimates of (6) remain finite-probe diagnostics only.

## 6. Frozen prospective finite probe

No shape science run preceded the registration. Exactly 8192 new shapes are drawn once, in eight blocks of 1024. Block b uses NumPy PCG64 with SeedSequence([202610277,b]), b=0,...,7, and 1024x64 Gamma(1/2,1) draws; P=Y/sum Y. One process is used, and BLAS/OpenMP thread limits are 1. No extra grid, new seed, adaptive batch, R draw, or automatic expansion is authorized. The source and proof are hashed before execution in registration.json.

At every shape the code records Z, all eight T_x, SPD eigen/Cholesky diagnostics, conditioning and factorization residual. Zero/nonfinite Gamma, nonfinite outputs, failed SPD or violated deterministic Z envelope halt unless the same fixed input is recomputed at 80 decimal digits. Condition number above 1e10 also invokes same-point 80-digit evaluation. No point is clipped, dropped, resampled, or reweighted. Wall time and CPU time are measured directly; all eight block summaries and the retained per-shape observables support the receipt.

The empirical full-moment cost uses U(g_hat,M_hat) with g_hat the minimum of **these eight sampled point estimates**, solely to test direction. It is neither a full-domain g estimate nor a C7 bound. The analytic M_Gamma endpoint is deliberately allowed to be coarse. The scientific outcome, source hashes and exact execution status belong to results.json, receipt.json and report.md.
