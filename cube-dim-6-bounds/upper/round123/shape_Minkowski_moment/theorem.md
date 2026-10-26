# R123: a complete boundary Minkowski budget from one determinant moment

Author proof version frozen before scientific execution, for primary review. Actual run status belongs to the separate report/receipt. Only the original body32 equal-weight Cauchy model is used; no new C6 bound is claimed and no lower-bound research is performed.

## 1. Normalization and the whole boundary

Let V consist of the32 vectors(1,s2,...,s6), sj in{−1,1}, with weights1. For independent S_i>0 satisfying E exp(−lambda S_i)=exp(−sqrt(lambda)), put Sigma=sum S_i v_i v_i^T and d=det(Sigma)^(-1/2). The exact representation, independently obtained by integrating the conditional Gaussian density with covariance2t^2Sigma, is

G(x)=C E[d Q_Sigma(x)^(-5/2)], Q_Sigma(x)=x^T Sigma^(-1)x,
C=Gamma(5/2)/(4pi^3)=3/(16pi^(5/2)).                 (1)

Indeed the conditional density is(4pi t^2)^(−3)d exp(−Q/(4t^2)). The substitution u=Q/(4t^2) gives integral_0^infinity t^(−6)exp(−Q/(4t^2))dt=16 Gamma(5/2)Q^(−5/2). Averaging its characteristic function exp(−t^2 sumS_i(v_i dot xi)^2) gives exp(−t sum|v_i dot xi|), so both the coefficient and the same Cauchy time are fixed without a factor2 ambiguity.

This is the same original Cauchy time and generator. Let q=6/5, F=[0,1]^5, I=integral_F G(1,y)^q dy. On the boundary of[-1,1]^6 use ordinary five-dimensional face measure dσ, with edges ignored as null sets. The original G is invariant under all coordinate permutations/sign changes because these permute the body lines with their equal weights. Consequently

integral_boundary G^q dσ =12*32*I=384 I.              (2)

A sampled Sigma need not have that symmetry. In particular no single-shape positive-face integral is multiplied by384 below.

For every SPD Sigma there is an exact full-boundary identity

integral_boundary Q_Sigma(x)^(-3) dσ(x)
  = integral_R6 exp(−Q_Sigma(z)) dz
  = pi^3 sqrt(det(Sigma)).                           (3)

For a direct proof, cubical ray coordinates z=rx on each face have Jacobian r^5 dr dσ. Nonnegative Tonelli and integral_0^infinity r^5 exp(−r^2 Q)dr=Q^(−3) give the first equality. The usual diagonalized Gaussian integral gives the second. This proof needs no single-shape coordinate symmetry.

## 2. Minkowski, its moment, and the full I upper bound

Applying Lq Minkowski on the entire boundary, (3) gives

||G||_Lq(boundary)
 <= C E[d ||Q_Sigma^(−5/2)||_Lq(boundary)]
 = C pi^(5/2) E[det(Sigma)^(−1/12)].                  (4)

Write M=E[det(Sigma)^(−1/12)]. Moment finiteness is proved below. To avoid assuming I finite while proving it, first restrict the expectation to the nested events E_n={n^(−1)Id<=Sigma<=n Id}. Apply Minkowski there, then let n increase. The positive integrands increase to(1); monotone convergence for their q-th powers gives(4) and finiteness of the full boundary norm. This also permits the known infinite G values on the trueS2, which form a face-null set. The estimate includes the entire t∈(0,infinity) potential. No short-time tail, middle-time term, far-time term, or CEH budget is omitted.

Using(2) only after the whole-boundary estimate,

I^(5/6) <= C (pi^3/384)^(5/6) M,
I <= I_M(M):= (3/16)^(6/5) M^(6/5)/384.              (5)

The pi factors cancel exactly. This is a complete upper bound on I, conditional only on a rigorous upper bound for the same model's M. An MC estimate of M cannot be substituted as a certified upper endpoint.

## 3. Optimality within the shape-weight Jensen family

Let w(Sigma)>0 be a scalar shape weight with finite positive expectation and the necessary mixed moment. Weighted Jensen followed by(3) gives

I <= C^q pi^3/384 * (Ew)^(q−1)
     * E[w^(1−q) det(Sigma)^(−1/10)].                 (6)

Hölder yields

(E det(Sigma)^(−1/12))^q
 <=(Ew)^(q−1) E[w^(1−q)det(Sigma)^(−1/10)].

Equality is obtained for w proportional to det(Sigma)^(−1/12); this recovers(5). Thus(5) minimizes this particular integrated Jensen expression, not the true I over all methods. The earlier determinant weight w=d yields the factor(E d)^(1/5). Since det(Sigma)^(−1/12)=d^(1/6), ordinary concave Jensen gives M^(6/5)<=(E d)^(1/5), so the new bound is at least as strong as that weight. It is strictly stronger when d is nonconstant almost surely; the actual positive iid scale law makes det(Sigma) nonconstant, for example by varying one scale on an open positive-density event. E d is finite: Sigma>=32 min_i S_i Id and negative moments of each S_i are finite. This comparison does not assert that either budget improves64.

## 4. Exact radial Rao–Blackwellization and a bounded shape observable

Put Y_i=1/(4S_i). These are iid Gamma(1/2,1). With R=sumY_i and P_i=Y_i/R, R~Gamma(16,1) is independent of P~Dirichlet(1/2,...,1/2). Set A(P)=sum v_i v_i^T/P_i and Z(P)=det(A(P))^(−1/12). Since Sigma=A(P)/(4R),

M = D E_P Z,
D=2 Gamma(33/2)/Gamma(16)
 =2*32!*sqrt(pi)/(4^16*16!*15!).                     (7)

There is a useful deterministic bound

0<Z(P)<=1/32,                                       (8)

with equality at the uniform shape. Here is a self-contained proof without enumerating basis subsets. Let c_J=det(V_J)^2 for six-column subsets J. Cauchy–Binet and V V^T=32Id give sum_J c_J=32^6. For α_J=c_J/32^6, differentiate det(sum_i t_i v_i v_i^T) at all t_i=1. Its derivative in t_i equals32^6 v_i^T(32Id)^(-1)v_i=32^6*6/32. Hence the α-weighted inclusion probability of each i is6/32=3/16. Weighted AM–GM now gives

detA =32^6 sum_J α_J product_{i in J}P_i^(−1)
 >=32^6 product_i P_i^(−3/16)
 >=32^12,

where product_i P_i<=32^(−32). At P_i=1/32, A=1024Id, so equality holds. This proves(8), finite first and second moments of Z, and M<=D/32<infinity independently of any p1 computation. It supplies no competitive certified expectation upper bound by itself.

Retaining the product before its maximum gives the additional accurate analytic endpoint

Z(P)<=32^(−1/2) product_i P_i^(1/64),
E product_i P_i^(1/64)
 =Gamma(16)/Gamma(33/2)*[Gamma(33/64)/Gamma(1/2)]^32,
M<=M_Gamma:=(2/sqrt(32))*[Gamma(33/64)/Gamma(1/2)]^32.  (8a)

The middle equality is the Dirichlet beta integral; its radial Gamma factor cancels exactly againstD. Thus substituting M_Gamma into(5) is already a rigorous complete analytic I bound. No decimal comparison of this particular AM–GM envelope with the competitive threshold is made here. The fresh probe observes the actual moment, and cannot be promoted to a certified replacement for M_Gamma without a separate upper-expectation proof.

For an independent analytic lower diagnostic, all diagonal entries of Sigma equal T=sumS_i, so Hadamard gives detSigma<=T^6. Since T has law1024S and E S^(−1/2)=2/sqrt(pi),

M>=1/(16sqrt(pi)).                                  (9)

The inequalities(8),(9) are analytic facts about this model. They are not C6 lower bounds.

## 5. Accurate g0 budget and general L1 interface

For a certified global g_L>0, the existing original-generator constant-truncation transfer gives the complete candidate operator bound

U_M(g_L,M_U)=6^(6/5) I_M(M_U)/g_L^(6/5)
            =(9/8)^(6/5) (M_U/g_L)^(6/5)/384.         (10)

This is a certificate expression, not an equality for the actual model cost. With the unproved target g0=1/65000, positivity of

S_M=64(1/(6*65000))^(6/5)
     −(3/16)^(6/5) M_U^(6/5)/384                    (11)

is exactly equivalent to the rational-power comparison

M_U^6 < (8/9)^6 *24576^5 /65000^6.                   (12)

Equivalently(E_P Z)_U < [(8/9)*24576^(5/6)/65000]/D, with D exact as in(7). Neither the g0 floor nor a rigorous expectation upper bound has been supplied. The existing121 global-g work remains necessary even if this complete I budget is promising.

For clarity about the L1 interface, homogeneity and cubical rays give, for c>0,

integral_R6 (G−c)_+ = (5/6)c^(−1/5)*384I.

Thus(5) pays integrability of the entire original positive-part kernel, including angular singularities and arbitrarily short/long times. If g_L is certified on the cube and0<c<g_L, division by64(g_L−c) gives the cube-majorizing normalized kernel; its mass is5 c^(−1/5) I/(g_L−c), minimized at c=g_L/6, yielding6^(6/5)I/g_L^(6/5). The accepted positive-source transfer for this same Cauchy generator then applies to general L1 inputs and all positive common cube radii. The new argument neither changes that generator/source proof nor mixes conditional Gaussian generators; it only replaces the complete I upper estimate and proves the needed L1 mass finite. No restriction to atomic inputs, one radius, or a truncated time window is introduced.

## 6. Frozen prospective experiment

The only proposed new science run observes Z(P), not109 point values or a finite-grid minimum. It uses8192 fresh iid shapes in8 fixed chunks of1024, one CPU with its own numerical threads limited to1, seed lists[20261027123,123,chunk_id]. It does not reuse109 shapes or seed streams. Conditional radial R is integrated exactly through(7).

The source, input SHA and this proof are submitted before execution. It will record Z means/empirical variance,8 block means, conditioning, exact formulas(7),(11),(12), and the empirical Minkowski budget relative to64. The latter only probes a direction: it is neither a deterministic expectation upper bound nor a rigorous budget exclusion. SPD/zero-Gamma/nonfinite or bound violations trigger same-input80-digit recomputation or a recorded stop, never clipping/dropping/reweighting/resampling. No automatic second run or sample expansion is authorized. A competitive rigorous upper bound for E Z is separate future work.
