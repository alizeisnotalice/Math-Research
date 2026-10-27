# C147: a finite positive two-dimensional pushforward for each C130 pair

This is a new upper-bound mechanism study in the unchanged C130 law. It does not replay the six-dimensional grid, evaluate a new k_h value or assert a new moment/C6 bound. The executed prototype is exactly one real cone/body pair; the subsequent m=3 and m=6 compiler packet is prepared only and requires a separate exact-hash release.

## 1. Preserve the joint law

In a fixed source cone C let E_1,...,E_6 be independent Exp(1). For a fixed unordered pair of original bodies read the exact nonnegative source rows b,c and set v_j=(b_j,c_j), Z=(A,B)=sum_j v_j E_j. No independence between A and B is assumed. Direct product integration gives, for s,t>=0,

    E exp(-s A-t B) = product_j (1+b_j s+c_j t)^(-1).

This derives from independent E_j, not independent A,B. Zero columns contribute exactly1 and may be removed. All other column lengths remain unchanged: normalizing their directions for fan bookkeeping is not permission to rescale the Exp variables. Column permutations preserve the law, and swapping both output coordinates preserves the k_h(A)k_h(B) integral. General direction-only matches are not law-equivalent.

The original C130 normalization remains

    mu2 = (46080/pi^3) sum_C alpha_C sum_(i<j) E[k_h(A_Cij) k_h(B_Cij)].

All alpha_C are the original positive rational weights. Thus replacing one six-dimensional expectation by its exact pushforward integral neither changes its coefficient nor changes the count moment, original fan or source law.

## 2. Rank2 initialization gives a positive density

Suppose the nonzero columns span R^2. Choose two independent columns u,w, reorder them counterclockwise so D=det(u,w)>0, and integrate the remaining variables later. The map (r,q)->u r+w q has Jacobian D. If rho is the unique rational vector satisfying rho.u=rho.w=1, the first two columns have density

    f_2(z) = D^(-1) exp(-rho.z) 1_{z in cone(u,w)}.

It is nonnegative, integrates to1, and has no atom on its boundary. Although a coordinate of rho can be negative, rho.z=r+q>=0 on this support. An expression must be used only on its designated sector. Every later independent ray variable changes the density by

    f_new(z) = integral_0^infinity exp(-r) f_old(z-v r) dr.          (1)

This is the convolution with the actual probability law of v E. Positivity and total mass1 follow inductively from (1) and Tonelli plus the translation z->z-vr. Absolute continuity follows because convolution of an absolutely continuous measure with any probability measure remains absolutely continuous. These arguments are the positivity and normalization certificate; no signed partial-fraction summand is declared positive.

## 3. Constructive finite exponential-polynomial closure

Maintain a finite angular fan of exact rational rays in the nonnegative quadrant. On each open sector C, write

    f_old(z) = sum_n P_n(z) exp(-lambda_n.z),

where P_n is a rational polynomial, lambda_n is rational, and the expression is allowed to have signed terms. It is only the complete sum that equals the positive density. Sector boundaries have planar measure0 and may use any consistent version; all expectation identities are unaffected.

For a term supported on C=cone(u,w), membership z-vr in C is equivalent to

    L_1(z)-r L_1(v)>=0,  L_2(z)-r L_2(v)>=0,
    L_1(z)=det(u,z), L_2(z)=det(z,w), r>=0.

For each L(v)>0 this supplies the upper candidate r<=L(z)/L(v); for L(v)<0 it supplies a lower candidate r>=L(z)/L(v). A zero denominator is not divided by: it instead requires L(z)>=0. Include0 among lower candidates. Thus the feasible r set is empty or a closed interval [ell(z),U(z)], whose bounds are maxima/minima of finitely many rational linear forms.

The interval is bounded for every nonzero v and z on the support: z-vr stays in the nonnegative quadrant, and any coordinate v_i>0 gives r<=z_i/v_i. One may include these redundant coordinate inequalities if desired. A zero v is handled separately and leaves f unchanged because integral exp(-r)dr=1.

The branch choices change only along the original fan rays and the new ray v. Indeed, equality between two nonzero candidates L_i(z)/L_i(v) and L_j(z)/L_j(v) defines a homogeneous linear form annihilating v; in two dimensions it is either identically0 or its zero line is span(v). Comparisons with the lower candidate0 and zero-denominator tests change on original boundary lines. The opposite halves of these lines lie outside the nonnegative support. Therefore adjoining v to the existing angular fan suffices. With d distinct nonzero column directions, there are at most d-1 nonempty interior sectors within the total support, hence at most m-1<=5 for m<=6 active columns. There is no extra fan arising from every endpoint-comparison pair.

For each new open sector choose one exact rational interior sample (the sum of its two normalized bounding rays). It fixes all comparisons; evaluate them there to identify the linear forms ell,U, or to identify the empty interval. These selections hold throughout that open sector, not merely at its sample.

Substitute a density term into (1). Put eta=lambda.v-1 and R(z,r)=P(z-vr). Its contribution is

    exp(-lambda.z) integral_ell(z)^U(z) R(z,r) exp(eta*r) dr.       (2)

R is a polynomial in r,z. If eta!=0, a polynomial antiderivative factor is

    Q(z,r) = sum_(k=0)^(deg_r R) (-1)^k eta^(-k-1) partial_r^k R,
    partial_r(exp(eta*r)Q) = exp(eta*r)R.

Evaluating at the linear ell,U gives two rational-polynomial times exponential-linear terms. If eta=0, integrate each monomial r^k as r^(k+1)/(k+1) and evaluate the same endpoints. Thus resonance is treated algebraically, without a singular limit or division by0. Terms with identical lambda in the same sector can be merged exactly, including cancellations.

This proves by induction a finite piecewise exponential-polynomial representation for every rank2 source pair, with polynomial degree at most m-2. Each convolution increases degree by at most1; the resonant branch is precisely how repeated rates/columns generate polynomial factors. Duplicate columns, parallel columns of different lengths, repeated directions and zero facet coefficients require no general-position assumption and no distinct-node partial fractions. The procedure is a constructive exact density generator, not a numerical approximation to a four-dimensional fiber integral.

## 4. An exact Laplace certificate for a generated density

For a sector with independent rational rays u,w and a term P(z)exp(-lambda.z), set z=u alpha+w beta. The Jacobian is |det(u,w)|. Expand P(u alpha+w beta)=sum p_ij alpha^i beta^j. For sufficiently large positive s,t, each term has convergent integral

    |det(u,w)| sum_(i,j) p_ij i! j!
      / [((lambda+(s,t)).u)^(i+1) ((lambda+(s,t)).w)^(j+1)].       (3)

Sum (3) over all terms and sectors, put over a common exact denominator, and verify that the numerator of the difference from

    product_j (1+b_j s+c_j t)^(-1)

is the zero polynomial. All operations are rational polynomial algebra. This checks the full density law, including resonant polynomial terms, rather than selected numerical density values. It is legitimate to verify first where all individual term integrals converge; the summed rational identity then determines the Laplace transform. At s=t=0 simplify the whole rational expression first, cancel removable factors and obtain mass1 from the product. Individual signed terms may have removable poles or may fail to be separately integrable there and must not be evaluated or clipped separately.

The constructive convolution proof supplies pointwise nonnegativity almost everywhere. The exact transform check independently detects algebraic generation errors; neither it nor a handful of float samples is a replacement for that proof. A future interval evaluator must enclose the full signed density expression before intersecting with the proved nonnegative range. Integrating positive parts of individual signed terms would change the measure and is forbidden.

## 5. Rank1 and rank0 are singular laws, not fake planar densities

If all nonzero columns are parallel, choose a fixed nonzero rational u in their direction and write v_j=d_j u with d_j>0. Then Z=u X with X=sum d_j E_j. Its scalar Laplace transform is product_j(1+d_j s)^(-1). Grouping lambda=1/d_j with multiplicities yields an exact partial-fraction representation in powers (s+lambda)^(-k), whose inverse is exp(-lambda x)x^(k-1)/(k-1)! times its rational coefficient. Equivalently use the one-dimensional convolution recurrence, with the same eta=0 polynomial branch. The full scalar density is nonnegative with mass1 by convolution; individual coefficients can be signed.

The correct expectation is the one-dimensional integral of k_h(u_1 x)k_h(u_2 x) against this scalar density. There is no planar Lebesgue density and no unmentioned arclength Jacobian. A zero coordinate of u is allowed. If all columns vanish, the law is the unit atom at(0,0). The actual Terra enumeration reports every current source pair as rank2, but this proof does not require that fact or silently omit the degenerate cases.

## 6. Executed one-real-pair prototype

Read primary C130 cone0, zero-based body pair(16,17). The exact rows are

    b=(0,0,0,0,1/60,1/30),  c=(0,0,0,0,1/60,0).

Four independent zero-column variables integrate out with total mass1. Let A=E5/60+E6/30, B=E5/60. The inverse is E5=60B, E6=30(A-B), Jacobian1800, giving

    f(A,B)=1800 exp(-30A-30B) 1_{A>=B>=0}.

The finite registered SymPy1.14.0 prototype actually verified the exact change of variables, positive support, normalization1, means E A=1/20 and E B=1/60, mixed moment1/900 and covariance1/3600>0. Thus A,B are explicitly dependent. It also verified the transform identity

    integral exp(-sA-tB)f dA dB
      =1800/[(s+30)(s+t+60)]
      =(1+s/30)^(-1)(1+(s+t)/60)^(-1).

The raw symbolic integral was returned by SymPy in a piecewise form; exact simplification of its difference from the rational product was checked to be0. The rational expression above is the clean normal form. Nonresonant integration-by-parts identities and the eta=0 monomial integrals were checked for degrees0..4, the maximum required at m<=6.

For nonnegative k_h, symmetry of exp(-30A-30B)k_h(A)k_h(B) on the full quadrant gives the additional special identity

    E[k_h(A)k_h(B)] =900 [integral_0^infinity exp(-30u)k_h(u)du]^2.

This follows by exchanging the two coordinates and splitting the quadrant into its two triangles. It does not assert independence. No k_h integral or special-function value was evaluated. Since0<=k_h<=p_h, truncating the triangular density to A<=R has omitted pair expectation at most p_h^2(2 exp(-30R)-exp(-60R)). This is a future analytic tail budget, not a completed quadrature certificate.

## 7. Optional analytic one-dimensional transform, with no value evaluation

For h>0 put p=erf(sqrt(h)) and I_lambda=integral_0^infinity exp(-lambda*a)k_h(a)da. The original C130 kernel formula is

    k_h(a)=[erfc(a/(2sqrt(h))-sqrt(h))
             -exp(2a)erfc(a/(2sqrt(h))+sqrt(h))]/2.

Its Gaussian derivative terms cancel exactly, giving

    k_h'(a)=-exp(2a)erfc(a/(2sqrt(h))+sqrt(h)),  k_h(0)=p.

Integration by parts gives lambda*I_lambda=p-J, where
J=integral exp(-(lambda-2)a)erfc(a/(2sqrt(h))+sqrt(h))da. A second integration by parts and completing the Gaussian square give

    (lambda-2)J=1-p-exp((lambda^2-2lambda)h)
                         *erfc((lambda-1)sqrt(h)).

Hence, for lambda!=0,2,

    I_lambda=[exp((lambda^2-2lambda)h)erfc((lambda-1)sqrt(h))
                -1+(lambda-1)p]/(lambda^2-2lambda).

All boundary terms vanish even for negative lambda: the probabilistic formula bounds k_h(a) by p*exp(a-a^2/(4h)), and erfc has the corresponding Gaussian tail. Differentiation and continuity in lambda are justified by the same integrable Gaussian envelope on compact parameter sets. The removable cases follow by differentiating numerator and denominator, or directly by the two integrations by parts:

    I_0=h(1+p)+sqrt(h)/sqrt(pi)*exp(-h)-p/2,
    I_2=h(1-p)-sqrt(h)/sqrt(pi)*exp(-h)+p/2.

In the executed pair only lambda30 is needed, so its unevaluated closed form is

    I_30=[exp(840h)erfc(29sqrt(h))-1+29p]/840,
    E[k_h(A)k_h(B)]=900 I_30^2.

This is an analytic identity, not a new numerical evaluation. Its signed numerator must be enclosed as a whole; no signed contribution is separately positive. No exp/erfc/erf value or incomplete-time integral is evaluated in the current prototype or the prepared compiler. A future finite numerical gate could evaluate this single expression, with strict whole-expression bounds.

## 8. Next finite packet and limits

The prepared next generator has exactly two new real patterns: cone0 pair(9,16), with3 active/3 distinct directions; and cone0 pair(0,1), with6 active columns, repeated columns and repeated directions. It constructs the fan and exact polynomial terms by (2), and checks (3) against the product after every generated stage. It evaluates no exponential values and no k_h integral. Its finite counts, expression caps, fixed cases and exact-hash release are separate from the executed m=2 prototype.

The positive-law reduction is general. Runtime, cancellation, number of polynomial terms, interval evaluation stability and a complete positive-integral improvement for all10416 source pairs are not yet certified. The9961 canonical source-pattern count is a bookkeeping fact, not an impossibility argument. Any later integration batch needs its own finite registration, exact density/source identity, full unbounded-tail accounting and independent audit. There is no new M upper, globalg or C6 result.
