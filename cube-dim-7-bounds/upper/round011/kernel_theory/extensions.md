# R011 theory-only extensions and sharp calibration

These are derivations, not additional numerical release. The default finite contract remains scalar DD tests, origin/high-symmetry values and256generic phases. In particular epsilon0 evaluation, fan moments and alternative quadrature do not start automatically.

## 1. A complete finite expression for G at a nonzero-node generic point

Suppose x is exact rational and b_r(Q^Tx)!=0 for every signed permutation Q and every catalog ray r. This is a finite exact integer test using the common-knot compiler. For each cone the time integrand product_i(t-i b_i)^(-1) is bounded near t0 and O(t^-7) at infinity. Its absolute time integral is finite. Hence the finite complete-cone sum can pass to epsilon0 termwise. Independently H_epsilon(x) increases to the true G(x), since the actual density p_t is nonnegative. Therefore

    G(x)=2c_x/(2pi)^7 sum_(y in O_x^+) sum_S alpha_S
            DD6[log|b|](b_r(y):r in S).                          (E1)

The order6 recurrence uses sorted exact rational nonzero knots. A repeated knot b has log|b| Hermite coefficient (-1)^(k-1)/(k b^k), k>=1. Although log|b| is not analytic through0, it is smooth near each repeated nonzero knot; recurrence convergence follows directly by epsilon-down-to0 and nonzero exact denominators. A mixture of positive and negative knots causes no problem for this limiting real recurrence. An equivalent complex implementation must keep the principal boundary argument at each sign, not replace log(-ib) by log|b| before taking the real part of a complex sum.

A zero node prevents this sufficient termwise argument. It does not imply that G is infinite: artificial cone singularities can cancel in the full expression. The exact singular set remains R001's S3. Record zero-node counts during the already released H pilot, but a new epsilon0 DD run needs a new release. A single generic G value is not a global face lower bound and is not an operator C7 counterexample.

## 2. Exact cone Fourier moments

Write xi=sum_i a_i r_i, a_i>=0, q_i=psi(r_i), alpha=|detR|/prod_i q_i and c_i=r_(i,j)/q_i. With h_m the complete homogeneous polynomial of degree m in seven variables,

    integral_cone xi_j^m exp(-epsilon psi(xi))/psi(xi) dxi
      =alpha m! h_m(c_1,...,c_7)
          /[(6+m)epsilon^(6+m)].                                (E2)

To prove this, use1/psi exp(-epsilon psi)=integral_epsilon^infinity exp(-t psi)dt, expand (sum_i a_i r_(i,j))^m by multinomial coefficients, and integrate each monomial against exp(-t sum_i q_i a_i). The factorial from each one-dimensional Gamma integral cancels its multinomial denominator; the result is alpha m! h_m(c)t^(-7-m). Integrating t gives (E2). For signed c and even m this remains nonnegative, as follows either from the left side or from m!h_m(c)=E(sum_i c_i E_i)^m for independent rate1 exponential E_i. On the original ordered cone all c_i are nonnegative.

Generate exact h_m by product_i(1-c_i z)^(-1). The full group signed-coordinate phase orbit has14elements and fibre46080. Thus the even coordinate moment, including the Fourier normalization, is

    J_(2k)=(2pi)^(-7) integral_R7 xi_j^(2k)
                   exp(-epsilon psi(xi))/psi(xi) dxi
      =720(2k)!/[pi^7(6+2k)epsilon^(6+2k)]
           *sum_S alpha_S sum_(j=1)^7 h_(2k)(r_(i,j)/q_i:i in S). (E3)

For k0 this gives840 A_sum/(pi^7 epsilon^6), exactly H_epsilon(0). Fourier differentiation gives (-1)^k partial_j^(2k)H_epsilon(0)=J_(2k). Signed permutations force the full second-moment matrix to J2 times the identity. For every real unit direction v,

    |D_v H_epsilon(x)|<=sqrt(J0 J2),
    |D_v^2 H_epsilon(x)|<=J2.                                   (E4)

The first uses Cauchy-Schwarz for the positive frequency measure, the second its isotropic second moment. Consequently the global Euclidean gradient norm is at most sqrt(J0 J2); there is no extra sqrt7 factor. These exact low-degree fan moments can calibrate future spatial certificates. They are not evaluated by this theory owner, and computing them is not permission to run a spatial cover or claim g.

## 3. Sharp Euclidean floor from the sealed spectral unit ball

Let P={xi:psi(xi)<=1}. The accepted R003 complete21spectral-vertex-orbit catalog is a different object from R009's46ordered-arrangement rays. A convex function attains its maximum on a compact polytope at some vertex: write any point as a convex combination and apply convexity. Hence the catalog plus exact comparisons

    max_vertex |xi|^2=1/2048

proves P lies in the Euclidean ball of radius1/sqrt2048. By homogeneity psi(xi)>=sqrt2048|xi|=32sqrt2|xi|. Equality holds along (1,1,0,...,0), since psi=64 and norm=sqrt2. This proves sharpness without invoking an external Khintchine theorem.

The primary thread actually performed the21exact Fraction comparisons as a new functional of the sealed catalog, not a fresh catalog enumeration. Its receipt is work/r011-primary-euclidean-floor.json, source catalog SHA256856e15aeaa785a72d4bbcd8d590141a8ad7b14fe179fae4671d5ef063f61e333. This owner read that receipt and derives the consequence, but did not execute the mathematical comparison. The source provenance is locked separately in sources.json. No repeat21-orbit run is needed merely to increase CPU usage.

For c=32sqrt2 and a real unit v, radial integration of the Fourier absolute bound yields

    ||D_v^k H_epsilon||_infinity
      <=Gamma(6+k)/(120pi^4 c^(7+k)epsilon^(6+k)).                (E5)

Indeed area(S6)=16pi^3/15, psi>=c|xi|, and the radial power is r^(5+k). The same domination proves holomorphy for |Im x|<c epsilon. The k1 improvement over c8 is (c/8)^8=2^20.

The sharper angular average for a coordinate square gives the optional calibration

    J0<=Gamma6/(120pi^4 c^7 epsilon^6),
    J2<=Gamma8/(840pi^4 c^9 epsilon^8).                          (E6)

The factor1/7 in J2 is the spherical mean of a coordinate square. If J0/J2 are later computed from (E3), compare their rigorous upper endpoints with (E6); a failure is a normalization or arithmetic diagnostic, never a reason to clip a kernel sum. Exact rational fan coefficients and the receipt-derived sharp floor have different provenance and should remain labelled.

## 4. Optional independent compact-time integral

The change t=1/u gives a finite-interval alternative to analytic DD integration:

    H_epsilon(x)=(2pi)^(-7) integral_0^(1/epsilon) u^5
           sum_fullcones alpha_S product_i(1-i u b_i)^(-1) du.  (E7)

The complete real integrand equals p_(1/u)(x)/u^2 and is nonnegative for real u>0. A single cone need not have positive real part. Formula(E7) uses the same accepted fan but an independent time-integration algorithm; one-cone quadrature cannot certify H.

For a real center u0 and complex disk radius rho with rho|x|/c<=1/2, every denominator in that disk has modulus at least1/2: at u0 it has modulus at least1, while the perturbation has size at most rho|b|<=rho|x|/c. Writing U=1/epsilon, a common holomorphic bound is

    M=(U+rho)^5 2^7 *645120*A_sum/(2pi)^7.

A Taylor polynomial of degree d on each real panel with halfwidth h<rho has total integrated remainder at most U M (h/rho)^(d+1)/(1-h/rho). This applies to the paired holomorphic rational expression, not to a nonholomorphic Re operation inside a complex Taylor implementation. One may take real coefficients after constructing the analytic pair.

If separately released, two inexpensive independent checks could use epsilon1/64 and degree64: e1 with rho4,h1,32panels; allones with rho1,h1/4,128panels. These safe radii even satisfy the older c8 floor. They are not part of the default R011 run. The bound above is a proposed deterministic error contract, not an achieved width or a promise that these160panels separate a target interval. Future release must lock code/parameters and test its actual enclosure; no extra quadrature starts from this document alone.
