# R126: the true first Gamma-count moment on the audited spectral cones

Pre-execution proof for primary review. This is a new auxiliary moment of the original body32 common-time Cauchy model, not a new C6 bound. No science has run at this proof stage.

## 1. Original law and positive Fourier representation

Let V be the32 unoriented body representatives(1,s2,...,s6), all weights1. Put Y_i iid Gamma(1/2,1), S_i=1/(4Y_i), Sigma=sum S_i v_i v_i^T, d=detSigma^(−1/2), psi(xi)=sum|v_i dot xi|, and N_h=sum_i1{Y_i<=h}. This count is under the independent, unnormalized Gamma law, not under Dirichlet coordinates and not a change of physical time.

For a>=0 define

kappa_h(a)=e^a E[e^(−a^2/(4Y))1{Y<=h}].

The full expectation is e^(−a), so0<=kappa_h<=1; kappa_h(0)=erf(sqrt(h)). The Gaussian integral identity pi^3 d=integral_R6 exp(−xi^T Sigma xi)dxi, independence, and positive Tonelli give exactly

J(h):=pi^3 E[d N_h]
 =integral_R6 e^(−psi(xi)) sum_(v in V)kappa_h(|v dot xi|)dxi. (1)

No Fourier inversion sign, oscillatory phase, or(2pi)^6 factor occurs in this Gaussian integral. Since sum kappa<=32, J(h)<=32R0, with the OLD exact R0=integral e^(−psi)=520066141/48769479475200. That R0 is read from002's producer and independent audit; it is not recomputed as a new peak. Equation(1) includes the whole original law and all small/large Y. It introduces no new physical Cauchy-time truncation.

## 2. One stable primitive, with derivative in h

For s>2, c=s−1>1, define L_h(s)=integral_0^infinity e^(−sa)kappa_h(a)da. Its integrand is nonnegative and bounded by e^(−sa). Substituting the Gamma density and swapping positive integrals gives

partial_h L_h(s)
 =e^(−h)/(sqrt(pi*h)) integral_0^infinity exp(−ca−a^2/(4h))da
 =T(c,h):=exp((c^2−1)h)erfc(c sqrt(h)).               (2)

This is the derivative with respect to h, not s. Completing the square proves(2). Moreover L_0(s)=0 by dominated convergence. Since

partial_h T=(c^2−1)T−c e^(−h)/(sqrt(pi*h)),
partial_h erf(sqrt(h))=e^(−h)/(sqrt(pi*h)),

the unique primitive is

L_h(s)=[T−1+c erf(sqrt(h))]/(c^2−1).                 (3)

It has initial value0; as h tends to infinity, T tends to0 and L tends to(c−1)/(c^2−1)=1/s, consistent with kappa_infinity=1.

Numerically use the stable form T=e^(−h)erfcx(c sqrt(h)). Local Arb has Tricomi U, so the implementation evaluates

T=e^(−h)U(1/2,1/2,c^2h)/sqrt(pi).                   (4)

To check the identity without an asymptotic approximation, U's positive integral gives U(1/2,1/2,x)=(2/sqrt(pi))integral_0^infinity e^(−xu^2)/(1+u^2)du. Writing(1+u^2)^(−1)=integral_0^infinity e^(−r(1+u^2))dr and swapping positive integrals yields e^x integral_x^infinity e^(−t)t^(−1/2)dt=sqrt(pi)e^x erfc(sqrt(x)). Thus(4) is exact and avoids a huge exponential times a tiny erfc.

## 3. Spectral nodes and their source

On one of the21 OLD fundamental cones write xi=sum_j r_j lambda_j, lambda_j>=0, psi(xi)=sum_j q_j lambda_j. Its mass is alpha_C=|det(r_1,...,r_6)|/product q_j. The integer fan input is pinned to SHA2567a57caad9b7713d9ce0896e51ab1dd402566feb1cba50829104c593a6d12a812. Its layout is the C016 parser's21*(det, six(q,r[6])) structure. The new parser binds each integer ray/rate/determinant by exact positive scaling to the already audited002 rational simplex and checks all q_j=psi(r_j), determinants, and weak signs. It does not rebuild the fan or run any old producer.

For each nonzero body label v, all v dot r_j must have one weak sign on that cone. This source property is explicitly checked in exact arithmetic. Only then

|v dot xi|=sum_j |v dot r_j|lambda_j
          =sum_j b_j E_j, b_j=|v dot r_j|/q_j,

under the normalized cone measure with independent E_j~Exp(1). These are spectral-frame nodes, not the signed physical-x evaluation nodes of an old Green kernel.

The tight identity V V^T=32Id gives psi(r)>=sqrt(32)||r||_2, while|v dot r|<=sqrt(6)||r||_2. Hence every positive b<=sqrt(3)/4<1/2. Thus c=1/b−1>1 and neither c=1 nor b=1/2 is approached as an actual pole. The exact input check enforces0<b<1/2. A zero b removes only its independent Exp variable, whose density integrates to1. If all b were zero, the remaining answer would be kappa_h(0)=erf(sqrt(h)), not1. In the actual cones m>=1 because the six rays span R6 and v is nonzero; this is also checked.

The SUM over all32 labels in(1) is hyperoctahedrally invariant. Therefore the full integral equals46080 times its fundamental-chamber integral. A single label term is not multiplied by that symmetry factor in isolation. We use46080*alpha_C, not C016's even-physical-kernel alpha2 convention. The exact existing normalization46080 sum_C alpha_C=R0 is a source-consistency check only.

## 4. Hypoexponential expectation and confluent jets

After deleting zero nodes, let b_1,...,b_m>0. For distinct nodes the density of sum b_j E_j, by elementary exponential convolution/partial fractions, is

sum_i b_i^(m−2) exp(−a/b_i)/product_(j!=i)(b_i−b_j), a>0.

Every term can be integrated against the bounded kappa_h, so

E kappa_h(sum b_j E_j) = [b_1,...,b_m]F_m,
F_m(b)=b^(m−2)L_h(1/b).                              (5)

This includes m1, where F_1=L_h(1/b)/b. For repeated nodes, the expectation is continuous by dominated convergence; F_m is analytic on0<b<1/2. Hence(5) extends to the ordinary Hermite confluent divided difference. Repeated nodes use F_m^(k)(b)/k!, never an artificial perturbation.

Set E=e^(−h)sqrt(h)/sqrt(pi), f=erf(sqrt(h)). Formula(3) simplifies to

F_m=A T+B f+C E+D,
A=b^m/(1−2b), B=b^(m−1)(1−b)/(1−2b), C=0,D=−A.

Since partial_b T=−b^(−2)(2chT−2E), one derivative sends

(A,B,C,D) to(A'−2chA/b^2, B', C'+2A/b^2, D').       (6)

All coefficients remain rational in b,h. Sympy differentiates only these rational functions, at orders0..m−1<=5. They are evaluated exactly as rational numbers at the exact nodes/h, then combined with Arb512 values of T,f,E and divided by k!. The Hermite table itself carries all signed subtraction at512 bits. No signed term is dropped.

Two exact algebraic consistency checks are registered: at h0 the jets vanish, and at h infinity F_m becomes b^(m−1). The latter's order(m−1) Hermite divided difference is exactly1 for every tuple, including repeated nodes. These are primitive/normalization checks, not old peak evaluations.

## 5. Finite exact moment certificate and comparison

Group identical positive-node multisets tau, with m=length(tau), and define the positive rational weights

w_tau=46080 sum_(C,v with tau(C,v)=tau) alpha_C.

There are at most21*32=672 ungrouped pairs, and sum w_tau=32R0. The exact observable is

J(h)=sum_tau w_tau [tau]F_length(tau).                (7)

The unique proposed run fixes h=1/4096,1/2048,...,1/32 and512 bits. Each functional gets an outward rational interval, then all positive weighted endpoints are summed with exact rational arithmetic. Dividing by an Arb enclosure of pi^3 yields an outward rational interval for the true mu_1(h)=E[dN_h]. Fixed384-bit outward grids are only endpoint representations; no precision increase, Monte Carlo, time grid, old point/native kernel or old peak producer is allowed. Any source, sign, special-function or precision failure is recorded and stops the run.

The old124 expression is B1*z, z=1−e^(−h), pi^3 B1=64R1, R1=67320683/3251298631680. Its complete saved lower/upper endpoints are copied from each124 h_i caps.b_lower[1]/beta_upper[1], with original-source hashes. These lower endpoints bound the OLD CAP EXPRESSION, not the true moment. A new upper endpoint below the old lower endpoint certifies a strict replacement-cap improvement. Fresh scalar enclosure of that same expression is checked to lie within the saved old interval; display decimals are never evidence.

The output changes only the k1 auxiliary moment cap. It does not optimize the125 joint-count model, pay higher count moments, prove a competitive M bound, certify globalg0, or change C6. Its use in the existing original-generator/alltime/generalL1 scheme leaves common Cauchy time and common cube-radius semantics unchanged.
