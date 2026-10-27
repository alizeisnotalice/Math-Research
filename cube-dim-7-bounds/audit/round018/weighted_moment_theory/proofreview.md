# R018 independent pure mathematical audit

2026-10-27. **Accepted within the explicit theoretical and finite-batch scope. No source erratum was found.** The determinant-weighted directional moment construction, two-node global lower function, exact root elimination and seven-constant extraction identities are correct. This audit does not certify that the lower function reaches a proposed floor, and does not promote the separate numerical batch to a globalg or C7bound. No mathematical numerical or symbolic program was run. Parent retains overall proof and integration responsibility.

## 1. Sources and dependencies

The entire frozen proof.md and registration.md were read, with their sources and receipt. The accepted R001 Gaussian/common-time identities, complete R010 fan, R016 origin normalization and original allr/generalL1 comparison are dependencies. Their declared bytes and the current theory manifest are hash-checked as metadata, not recomputed as mathematical certificates. Source counts such as212top simplices and46rays are retained dependencies. The independent clone HEAD is5b24e19fc244cea2589469206bbfcb394f48ea9d, while the source theory itself records its earlier base14f3507f406185667f84146328f1db9da3f9949a. No assumption that these are the same checkout is made.

Only audit/round018/weighted_moment_theory was written. Source documents, earlier sealed artifacts and other workers' files were not changed. The separate seven-exact/three-scalar execution belongs to its numerical owner; this audit did not read a later numerical result as a substitute for proof, execute its program or extract the coefficients itself.

## 2. Time-one density, tilt and directional moments

Let P~Dirichlet(1/2)^64 and R~Gamma(32,rate1) be independent, A=sum_v vv^T/P_v, B=A^-1 and d=(detA)^(-1/2). The original stable mixture has Sigma=A/(4R). Conditional on R,P, the time-one Gaussian density is

    pi^(-7/2) R^(7/2)d exp(-R Q_x), Q_x=x^T Bx.

Positive R integration gives

    p1(x)=K0 E_P[d(1+Q_x)^(-a)],
    a=71/2, K0=Gamma(a)/(pi^(7/2)Gamma32).

The full-time potential obtained from the same original common-time mixture is

    G(x)=KD E_P[d Q_x^(-3)],
    KD=Gamma(a-3)/(pi^(7/2)Gamma32).

Here p1 is the actual time-one density. At0 the positive Fourier/Gaussian integral identifies its value directly, avoiding a selected-point issue from merely a.e. equality of density versions. It is finite, positive, and p0=K0 E d. G(0) is instead infinite. These two objects must not be interchanged.

Since d is finite and positive on the simplex interior, the normalized law dnu=d dPr/(E d) is a genuine probability with the same full interior support as the Dirichlet law. Gamma recurrence gives KD/K0=1/b with b=(65/2)(67/2)(69/2), hence

    G(x)=(p0/b) E_nu Q_x^(-3).

The tilt is by d alone after radial integration, not by a raw radial moment or the trace proxy. The stated normalization is correct.

For the original symbol psi, define the **density** Fourier moment

    J_(2k)(x)=(2pi)^(-7)integral (x dot xi)^(2k)exp(-psi(xi))dxi.

These moments differ from the R011 H moments, whose integrand has an extra1/psi. The positive psi floor supplies integrability of all required polynomial frequency moments. Also A>=2048I gives Q_x<=|x|^2/2048 uniformly. Differentiating

    p1(tx)=p0 E_nu(1+t^2Q_x)^(-a)

near t0, justified by this uniform Q bound, yields

    mu_k=E_nu Q_x^k=k! J_(2k)/[(a)_k p0 (2k)!].

Thus the three divisors are2a,12a(a+1),120a(a+1)(a+2). The factors are density-Taylor factors; no potential time-integral divisor may be imported. This independently verifies the central normalization.

## 3. Nondegeneracy for every nonzero x

For x!=0, Q_x>0 at every interior P, and all its positive moments are bounded. At uniform P, the body sign cancellation gives sum_v vv^T=64I, so A=4096I and Q_x=|x|^2/4096>0. Along the interior path P_i=1-63delta, P_j=delta for j!=i, the divergent part of A is delta^-1(64I-v_i v_i^T). Since |v_i|^2=7, this Gram matrix is SPD with least eigenvalue57. Therefore Q_x tends to0 for every nonzero fixed direction x.

Continuity of Q on the interior and full support of nu mean every interior path value lies in its scalar support: a neighborhood of the shape maps into a neighborhood of that value and has positive nu measure. The positive values approaching0 cannot form a two-point support. This avoids incorrectly assigning positive probability to the individual uniform shape itself. Consequently, with m=mu1,

    V=mu2-m^2>0,
    W=E[Q(Q-m)^2]=mu3-2m mu2+m^3>0,
    H=m mu3-mu2^2>0.

The last inequality is strict Cauchy-Schwarz for sqrt(Q) and Q^(3/2); equality would force Q constant on its positive support. This proof applies to every nonzero x, including e1 and the cube corner even when trueG is infinite there. No nonseparation in floating intervals can be interpreted as model degeneracy.

## 4. Positive two-node rule and global Hermite lower

Set c=W/V, S=m+c and T=H/V=mc-V. The identity H=mW-V^2 follows by expansion. The polynomial z^2-Sz+T is orthogonal to1,z; this follows by substituting the three moment definitions. Its discriminant is(m-c)^2+4V>0, its sum of roots is positive and product T>0, so both roots r1<r2 are positive. Its value at m is-V<0, hence r1<m<r2 and the stated weights are strictly positive and sum1. The orthogonal-polynomial sign-change argument also places them within the convex support interval. Cubic division by this quadratic proves exactness for all degrees0through3.

For f(z)=z^-3, f''''(z)=360z^-7>0 on the whole positive half-line. The cubic interpolating f and f' at both nodes has pointwise remainder

    f''''(theta)(z-r1)^2(z-r2)^2/4! >=0

for every positive z, because the interval containing z and the nodes stays positive. It is a global lower Hermite polynomial, not a local expansion near the nodes. The polynomial expectation is finite since the first three moments exist. The inverse moment may be infinite; integrating a pointwise nonnegative error above an integrable polynomial is still legitimate. Thus

    E_nu Q^-3>=lambda1/r1^3+lambda2/r2^3.

The inequality is strict for this infinitely supported law, as there is positive mass outside the two interpolation nodes. No quantitative positive remainder is produced or required. For other laws, degeneracies need their own treatment; the model's strict V,H proof is what licenses this formula here.

## 5. Root-free inverse-third expression

Let the Jacobi matrix be J=[[m,sqrt(V)],[sqrt(V),c]], with determinant T>0. The rule's nodes and weights are its eigenvalues and first-coordinate spectral weights. Consequently its inverse-third moment is the(1,1)entry of J^-3. Direct algebra using J^-1=T^-1[[c,-sqrt(V)],[-sqrt(V),m]] gives

    gauss_inverse3=[c^3+(2c+m)V]/T^3
                 =[W^3+2WV^3+mV^4]/H^3.

Both equalities and every power of V are correct. In particular the two terms after W^3 are2WV^3 and mV^4. All numerator terms and H are strictly positive in the actual model. Jensen for the positive two-node probability measure also gives gauss_inverse3>=m^-3.

Thus L2=(p0/b)gauss_inverse3 is a valid positive lower function for G. For exact rational directional moment coefficients and rational x, all of m,V,W,H and gauss_inverse3 are rational; pi^-7 appears only through p0. Roots, near-root subtraction and quadrature-weight interval division are unnecessary. The partition coefficient R_part2 and this inverse-third function must remain distinct symbols, as the source already requires.

## 6. Six mixed Fourier coefficients and all combinatorial factors

B7symmetry of exp(-psi) removes odd coordinate powers and makes each mixed even moment depend only on its partition. The tilted law has the same symmetry as well: the body-line permutation conjugates A,B and preserves determinant weight d. For k<=3 the partitions are(1),(2),(11),(3),(21),(111). Expanding(x dot xi)^(2k) gives multinomial coefficients6for the two-square fourth-order term,15for the ordered fourth/two sixth-order term, and90for three squares. Therefore the displayed J2,J4,J6 formulas, including their p1,p2,p3 versions, are correct. No extra factorial from the Taylor normalization belongs in these Fourier monomial coefficients.

On an ordered fan cone put u_r=xi_r/q_r and alpha=|detW|/product(q). Positive exponential integration reduces polynomial moments to seven independent rate-one exponential coordinates E_r:

    ell_j=sum_r(u_r)_j E_r,
    E product_r E_r^(n_r)=product_r n_r!.

This is exact monomial integration, not sampling. For a monomial on l fixed even-powered coordinates, the group sum has factor2^7(7-l)! times the sum over ordered distinct coordinate assignments. Fourier normalization divides out2^7 and leaves(7-l)!/pi^7. It follows that

    R0=7! sum alpha,
    R1=6! sum alpha E P2,
    R_part2=6! sum alpha E P4,
    R11=5! sum alpha E(P2^2-P4),
    R3=6! sum alpha E P6,
    R21=5! sum alpha E(P2P4-P6),
    R111=4! sum alpha E(P2^3-3P2P4+2P6).

Here R_lambda=pi^7T_lambda, R0=pi^7p0, and P_(2j)=sum_i ell_i^(2j). In the repeated-index cases the differences are already **ordered** distinct sums: P2^2-P4 includes both orders of a two-index pair, and the cubic combination includes all six orders of a distinct triple. An additional division by2or6 would be erroneous. For partition21 the two orders of assigning unequal exponents are already included in P2P4-P6. The source formulas have exactly the correct5! and4! factors. Also R0=5040A_sum agrees with R016.

All seven constants are positive rational cone integrals; rho=R_lambda/R0 is rational. With common L=lcm(q), each ell has one factor L^-1, so a homogeneous degree2k contraction divides by L^(2k) exactly once. Factorial contraction does not change this degree. The source sparse schedule has the stated homogeneous support bounds and finite pair-product expression: seven linear/square/fourth/sixth chains and the P2^2,P2^3,P2P4 chains. The count is an a priori algebraic schedule bound, not executed operations or timing evidence in this audit. No full group unfolding, higher moment degree or numerical integration is necessary for these identities.

## 7. Allcube, essential-floor and point-gate quantifiers

For every x!=0, the polynomial moments are smooth, V,H are strictly positive and the displayed rational lower function is smooth and positive. Its B7invariance follows from the moment formulas. Scaling x by t multiplies m,V,W,H by t^2,t^4,t^6,t^8 respectively, so gauss_inverse3 and L2 scale by t^-6. At0 the chosen extended value is infinity; substituting0/0 into the rational expression is not allowed.

The boundary of the cube is compact and excludes0, so L2 has a positive minimum there. For any nonzero interior point x=r y with y on the cube boundary and0<r<=1, L2(x)=r^-6L2(y)>=L2(y). Thus proving g<=inf_boundary L2 suffices on the entire cube, not merely the face. Symmetry reduces equivalent faces, but no conditional shape symmetry is assumed before averaging.

Conversely a fixed-point upper enclosure for L2 below a target is only an obstruction to reaching that target **using this lower function**. It is never an upper enclosure for actualG. Continuity gives failure on a relative-open positive-face-measure neighborhood even if the failing point is a face corner; this also rules out an essential-face certificate from this function. The original G at a corner or e1 may be infinite. Therefore such a gate does not refute the true Green floor. Passing the three predeclared points is only a necessary diagnostic, not a global certificate.

On x1=1, the six variables s_i=x_i^2 in[0,1] enter through three power sums u1,u2,u3. The two stated moment matrices are necessary PSD constraints, obtained by summing s(a+bs)^2 and(1-s)(a+bs)^2. Their closed feasible set is bounded: its diagonal and determinant inequalities bound0<=u3<=u2<=u1<=6, so this is a compact outer set. It allows measures beyond six equal-weight atoms, so failure there need not be failure on the actual domain. For a target g, multiplying the lower inequality by the actual positive bH^3 gives precisely

    p0(W^3+2WV^3+mV^4)-b g H^3>=0.

Coverage, sign conditions and any later polynomial/interval/SOS certificate must be proved separately. The full-support proof guarantees H>0 on actual nonzero physical directions, not automatically on every artificial outer-set triple. No optimization or global proof follows merely from the formal three-variable reduction.

## 8. Strict scope and stop

Accepted: the pointwise density normalization, determinant tilt, exact density Fourier moments, full-support nondegeneracy, global Hermite lower including infinite inverse moments, root-free formula, six mixed coefficient factors, seven finite extraction identities, denominator scaling and the stated allcube/point-gate implications. No source correction is needed; source proof and registration hashes accompany this report.

Not claimed by this audit: computed moment coefficients, success or failure of a numerical gate, certified globalg, C7, numerical lower remainders, source execution timing, or feasibility of future global polynomial work. The ten-output convention is seven exact coefficients plus three pi-scaled lower-function values, not ten transcendental library calls. The separate numerical owner must preserve the frozen inputs, exact signs and the R016 outward target provenance, and must not interpret a surrogate upper endpoint as a Green upper bound.

No mathematical numerical or symbolic software, samples, phase/DD, MATLAB, commit or push was used. Only full proof/spec reading, independent mathematical reasoning and read-only metadata hashes were performed. This pure audit is complete and stops; no follow-on is started.
