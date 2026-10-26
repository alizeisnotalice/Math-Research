# A057: positive-mixture log-curvature and a whole-box interface

All statements below concern the original unit-weight body32 model only:

    V={(1,s2,...,s6): si in {-1,1}},  A=[v:v in V],
    f_t(z)=t/[pi(t^2+z^2)],  p_t=density of sum_v v Z_(v,t),
    Z_(v,t) independent with density f_t.

Fix either I=[a,T] with 0<a<T<infinity or I=[a,infinity) with a>0. Set F(x)=integral_I p_t(x)dt. Thus F is exactly B_(a,T) or H_a. This is not full G; the lower cutoff is indispensable. All bounds hold on all R6, without chamber, projection-separation, or distance-to-S3 hypotheses.

## 1. Six-factor representation, finiteness, and differentiation

Choose any six distinct body columns forming an invertible matrix N and write A=[N,M], reordering independent variables if needed. Let nu be the product law of 26 standard Cauchy variables z. A change of variables in the six selected variables gives the everywhere-defined continuous density

    p_t(x)= integral k_(t,z)(x) nu(dz),
    k_(t,z)(x)=|det N|^-1 product_(i=1)^6 f_t((N^-1(x-tMz))_i).

This is ordinary positive convolution, not conditioning on a zero-measure event. In particular

    0 < p_t(x) <= |det N|^-1 pi^-6 t^-6,
    0 < F(x) <= |det N|^-1 pi^-6/(5a^5).                    (1)

The same upper bound covers the finite window. It pays the entire infinite time tail without a numerical cutoff or a missing remainder.

For l_t=log f_t and real z,

    l'_t(z)=-2z/(t^2+z^2),
    l''_t(z)=2(z^2-t^2)/(t^2+z^2)^2,
    |l'_t|<=1/t,       -2/t^2<=l''_t<=1/(4t^2).           (2)

The nonnegative gap identities are

    1/t^2-(l'_t)^2=(z^2-t^2)^2/[t^2(t^2+z^2)^2],
    l''_t+2/t^2=2z^2(z^2+3t^2)/[t^2(t^2+z^2)^2],
    1/(4t^2)-l''_t=(z^2-3t^2)^2/[4t^2(t^2+z^2)^2].

For b=N^-1 u, define sigma=D_u log k and tau=D_u^2 log k. Then

    |sigma|<=||b||_1/t,
    -2||b||_2^2/t^2<=tau<=||b||_2^2/(4t^2),
    D_u k=k sigma,     D_u^2 k=k(sigma^2+tau).             (3)

The absolute first and second derivatives are at most the peak in (1) times, respectively, ||b||_1/t and (||b||_1^2+2||b||_2^2)/t^2. These uniform bounds are integrable against dt nu(dz) on either I because a>0. Mixed derivatives follow directly or by polarization. Dominated differentiation gives F in C^2 everywhere (and higher orders can be treated likewise). No full-G regularity is being inferred.

Normalize the positive measure k_(t,z)(x)dt nu(dz) by its finite strictly positive total F(x). Expectations under this probability law are denoted E_x. Direct quotient differentiation gives

    D_u log F=E_x sigma,
    D_u^2 log F=E_x tau+Var_x(sigma).                      (4)

Consequently, using t>=a and ||b||_1^2<=6||b||_2^2,

    -2||N^-1u||_2^2/a^2 <= D_u^2 log F(x)
       <= (25/4)||N^-1u||_2^2/a^2.                        (5)

The upper constant is valid but is not asserted optimal. The lower inequality, not the upper one, pays the lower tangent certificate. In matrix form its left side is the fixed Loewner bound Hess log F >= -2 N^-T N^-1/a^2.

## 2. A stronger fixed-fiber version for the same model

This improvement is proved independently; it does not transfer data between models. Let R be any fixed 32-by-6 right inverse AR=I. Partition R=(R_N;R_M) according to [N,M]. In the positive full-rank formula

    p_t(x)=|det N|^-1 integral_(R26)
            [product_i f_t((N^-1(x-My))_i)] [product_j f_t(y_j)]dy,

make the x-dependent translation y=R_M x+w. Since NR_N+MR_M=I, this becomes an integral against the fixed Lebesgue measure dw:

    p_t(x)=|det N|^-1 integral product_(i=1)^32 f_t((Rx+Kw)_i)dw,
    K=(-N^-1M; I26).                                      (6)

The translation Jacobian is one. Thus (6) is a positive mixture on the fixed space I times R26 even though the unnormalized measure there has infinite total mass. Its weighted total F(x) is finite by (1), which is all normalization requires.

Derivative domination in this representation needs a separate argument. Integrating |l'_t|<=1/t along a real shift gives

    product_i f_t((R(x0+h)+Kw)_i)
      <= exp(||Rh||_1/a) product_i f_t((Rx0+Kw)_i).

On any bounded neighborhood of x0, this factor is uniformly bounded. First and second derivative ratios are bounded by ||Ru||_1/a and (||Ru||_1^2+2||Ru||_2^2)/a^2. The right-hand product is integrable in (t,w), since its integral equals |det N|F(x0). Dominated differentiation is therefore legitimate locally at every x0, including for the infinite tail. Repeating (4) with 32 factors proves

    Hess log F(x) >= -2 R^T R/a^2.                        (7)

For body32 the sign cancellations give AA^T=32I6. Use the single, direction-independent right inverse R0=A^T/32. Then

    R0^T R0=I6/32,
    D_u^2 log F(x) >= -||u||_2^2/(16a^2),
    log F(c+u) >= log F(c)+grad log F(c).u
                              -||u||_2^2/(32a^2).         (8)

This is stronger than the specified six-factor bound and uses the same F and gradient data. It has no positive all-time moment J2 and does not contradict A040's S3 obstruction: a stays strictly positive, and the loss grows as a^-2.

Within the class of bounds (7) obtained by a fixed right inverse, R0 is Loewner-optimal. Write any R=R0+D; then AD=0, R0^TD=0, and R^TR=I/32+D^TD >= I/32. This does not claim the actual Hessian of log F attains the bound or that no other method can improve it.

## 3. A concrete legal six-body basis and its face quadratic

Let J6 be the all-ones matrix, D=diag(-1,1,1,1,1,1), and

    N=(J6-2I6)D.

Its first column is (1,-1,-1,-1,-1,-1), and column j>=2 has first entry 1, entry j equal to -1 and all other entries 1. Thus these are six distinct columns of body32. Since J6-2I6 has eigenvalues 4,-2,-2,-2,-2,-2, det N=128. Moreover

    NN^T=4I6+2J6,
    N^-T N^-1=(NN^T)^-1=I6/4-J6/32.

On x1=1 write u=(0,d), d in R5. The exact face quadratic is

    Q6=I5/4-J5/32,
    q6(d)=d^T Q6 d=(1/4)sum_i d_i^2-(1/32)(sum_i d_i)^2. (9)

Its eigenvalues are 3/32 in the constant direction and 1/4 in the four orthogonal directions. For the full-frame bound (8), Q32=I5/32. The difference Q6-Q32=(7I5-J5)/32 is positive definite. No optimality over six-column bases is claimed.

For a centered face box |d_i|<=h_i, the exact worst q6 loss is

    max q6=(1/4)sum h_i^2-(1/32)min_(eps_i=+-1)(sum eps_i h_i)^2.

For equal halfwidth h, max q6=39h^2/32, whereas max q32=5h^2/32. The coefficient ratio is exactly 39/5. This is only a curvature-cost comparison, not a measured coverage speedup.

## 4. The exact whole-box tangent bound

For either fixed Q in (9) or Q32, let c=(1,c2,...,c6), and let the target box have free-coordinate intervals [l_i,r_i]. It may cross chamber boundaries and projection-collision sets. Put v=grad_face F(c)/F(c), and let Dbox be the displacement rectangle product_i[l_i-c_(i+1),r_i-c_(i+1)]. Integrating the lower Hessian along the whole segment gives

    F(c+(0,d)) >= F(c) exp(v.d-d^TQd/a^2).                (10)

The factor 1/2 from the Taylor integral is already included in Q/a^2. Because Q is positive semidefinite, v.d-d^TQd/a^2 is concave. A concave function on a rectangle is bounded below by its smallest vertex value, using the rectangle's convex-combination coordinates. Thus, with Vbox its 32 corners (duplicates harmless),

    inf_box F >= F(c) exp(min_(d in Vbox)
                                      [v.d-d^TQd/a^2]).  (11)

This retains the joint directional-gradient and quadratic cost at each corner. Subtracting separately max quadratic and max absolute gradient is valid but can lose margin. The theorem holds for the whole box, hence also its intersection with the closed ordered face or A040's paid compact domain. It does not itself prove a collection of boxes covers either domain.

Only one center value and five face derivatives are needed. Twenty-one point Hessian fields from A048 are unnecessary here; conversely, those fields by themselves do not supply a whole-box remainder.

### Rigorous interval inputs, retaining the shared value denominator

Suppose the same-window center jet has certified rational endpoints

    0<L<=F(c)<=U,      g_i^-<=partial_(i+1)F(c)<=g_i^+.

At each corner d set

    b_d=sum_i min(g_i^- d_i,g_i^+ d_i),
    f_d=min(U,max(L,b_d)),
    r_d=b_d/f_d-d^TQd/a^2,
    C_d=f_d exp(r_d).                                    (12)

Then inf_box F >= min_d C_d. Indeed, for fixed positive f the gradient interval minimizes its dot product at b_d. The derivative of log f+b_d/f is (f-b_d)/f^2, so its exact minimum over [L,U] occurs at the displayed clamp f_d, including b_d<=0. Minimizing over the interval jet and the rectangle commutes with the finite corner minimum. Formula (12) is the exact robust minimum of the tangent bound over the independent endpoint box; any true correlations can only improve it.

The executable checker computes rational b_d, f_d, q_d and r_d exactly, evaluates the exponential with Arb directed balls, and accepts only if every outward lower endpoint is >= the rational target. An inconclusive endpoint returns UNRESOLVED. Matching model, center, and time-window metadata are checked structurally; the checker cannot manufacture or independently authenticate the supplied kernel-jet proof. Its PASS is explicitly conditional on that proof.

## 5. Multiple bases, scales, singular strata, and coverage scope

All bounds apply for every a>0 and every T>a with the same unit-weight model. There is no fixed-a spatial homogeneity assertion: the exact scaling is

    H_a(r x)=r^-5 H_(a/r)(x),
    B_(a,T)(r x)=r^-5 B_(a/r,T/r)(x),  r>0.

Therefore changing a or T requires matching new or rigorously transformed value/gradient data. Full-G jets, a different time mixture, sparse4/augmented-frame data, or a band value with an H gradient cannot be substituted. After a completely paid face certificate, full G>=F and the genuine full-G signed-permutation symmetry/homogeneity may transfer the floor to the cube, as in A040. That future coverage step has not been done here.

Zero and repeated Fourier projections are irrelevant to this real-space proof. It neither divides by their differences nor excludes them. A future center evaluator must nevertheless evaluate their legitimate Hermite limits, as in A048. S3, S2, ties and zeros all lie in the domain of this smooth cutoff theorem; the uncut Green function may behave differently.

For multiple six-body bases, each entire parabola is valid. One may certify a complete box with each fixed Q and take the largest certified whole-box lower bound. Do not take the pointwise minimum of several quadratics and apply one 32-corner concavity rule: that minimum need not be convex. For example q_+=(d1+d2)^2 and q_-=(d1-d2)^2 have min(q_+,q_-)=0 at all four corners of [-1,1]^2, but value 1 at (1,0). Thus the corresponding negative minimum has a lower interior/edge value than all corners. Adding a small common positive quadratic preserves the counterexample. The full-frame improvement (8) avoids this problem because it uses one fixed matrix, proved before optimizing any direction.

## 6. Acceptance and the next permitted experiment

This round proves the requested six-body lower bound, its optional 25/4 upper bound, and the stronger fixed-right-inverse lower bound, with an executable whole-box interface. Seven fixed symbolic/matrix/interface configurations are the complete batch. Synthetic interface jets are not body32 observations.

A later pilot would need new matching H_a or B_(a,T) values and five face gradients at registered centers covering ordinary low-value regions, S3 strata away from the already paid S2 tube, tie/zero boundaries, and seams with the tube. At least one difficult S3 anchor can be chosen from A040's family (1,1,1,1,b,c), with exact rational b,c and a>0 chosen before execution. Corners may cross all projection strata because this theorem pays them analytically. No such jet or kernel call is run here; a finite point pilot still would not be a global cover.
