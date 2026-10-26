# A051 all-middle independent proof review

We fix the original body32 common-time Cauchy model, rho=1/64 and epsilon=1/1000. A contains coordinate1, 1≤k=|A|≤5, B=Aᶜ, m=6−k, v=ones, w=1_A−1_B, U_A=span(v,w). The domain is exactly

    D_A={x=(1,y), y∈[0,1]^5 : dist(x,U_A)<rho,
                                  mean_B(x)∈[1/4,3/4]}.

We accept the previously reviewed033 projection peaks,043 scalar Cauchy-average/complete-triple certificate, and049 joint-marginal theorem as frozen inputs. This review checks their applicability to all five types; it does not rerun their writers, old budgets or enumerations.

## 1. Actual perturbation and all signed nearest lines

Write x_A=alpha1_A+e, x_B=u1_B+z, with e,z zero-sum and e1=1−alpha. For x0=1_A+u1_B,

    h=x−x0=(e−e1 1_A,z),
    |h|²=|e|²+|z|²+k e1².

The coordinate functional on the k-block zero-sum space has squared norm(k−1)/k. Thus k e1²≤(k−1)|e|² and |h|≤sqrt(k)dist(x,U_A)<eta_k=sqrt(k)rho. For k=1 the A residual is zero and the bound remains valid. No factor sqrt(3) is reused indiscriminately across types.

For a vector with positive value p on A and value q on B, every signed body scalar product has absolute value≤kp+m|q|. This value is attained by v if q≥0 or by w if q<0. Since every body column has norm sqrt(6), the exact squared distance to the entire union of body lines is

    kp²+mq²−(kp+m|q|)²/6=(km/6)(p−|q|)².                 (1)

This is an all-signed-lines statement. For |a_v|<1/2, p=1−a_v≥1/2 and q=u−a_v. If q≥0 then p−|q|=1−u≥1/4. If q<0 then p−|q|=1+u−2a_v≥u≥1/4. For |a_w|<1/16, p=1−a_w≥15/16, q=u+a_w≥3/16 and p−q=1−u−2a_w≥1/8. Norms are at least sqrt(k)p. The distance-to-set and norm functions are1-Lipschitz, so the following bounds hold for actual x after subtracting eta_k:

    hv0=sqrt(k)(1/2−rho),
    hv1=(1/4)sqrt(km/6)−sqrt(k)rho,
    hw0=sqrt(k)(15/16−rho),
    hw1=(1/8)sqrt(km/6)−sqrt(k)rho.                         (2)

All are strictly positive. The nontrivial weakest reduced gates are m/384−rho²>0 and m/96−rho²>0; m≥1 suffices at rho=1/64. The other two gates are rational. This proves positivity without changing thresholds or repairing a type.

## 2. All other signed pair planes at original x

At x0, a different unsigned two-block partition necessarily has a block mixing at least one coordinate of value1 with one of valueu. Its variance is at least(1−u)²/2, because minimizing just those two squared deviations already gives that amount. Thus the plane distance is≥(1−u)/sqrt(2)≥1/(4sqrt(2)). A signed partition has no larger projection gain than its unsigned version at x0≥0, so the same lower bound applies.

If the unsigned partition is A/B itself, a different signed plane must have opposite signs within some block of size≥2. For such a pair, (e_i+e_j)/sqrt(2) is perpendicular to that signed block line and its inner product with x0 is sqrt(2) times that block's value. The distance is therefore≥sqrt(2)u≥sqrt(2)/4. A singleton block causes no exception: mixing signs can only occur in the other block. Constant signs in both blocks give U_A itself.

Subtracting the actual perturbation bound gives, for every body pair plane U≠U_A,

    dist(x,U)≥g_k=1/(4sqrt(2))−sqrt(k)rho>0.                (3)

The gate1/32−k rho²>0 holds for all k≤5. The only body lines in U_A are v,w: their block values must be±1 and the first coordinate fixes the A sign. Hence any pair avoiding v, or avoiding w, cannot span U_A, so every pair used below is covered by(3).

Crucially (3) is a bound for original x. It is not asserted for x−a b uniformly over a. Third-distance zeros along that line are retained and integrated.

## 3. True pair coarea and the033 R2 type

The actual selected-column Gram determinant is36−(k−m)²=4km, and J2=2sqrt(km). On the perpendicular quotient, deleting v,w leaves the original restricted symbol unchanged. The033 normalized table states

    Q4,1(0)/J2=R2,min(k,m)/(2π)^4.

The type uses min(k,m) because a quotient plane is unchanged by changing a spanning vector's global sign and is invariant under full signed permutations. This symmetry does not exchange the positive-face roles of A and B in(2) or in the eventual counts.

The two large actual parameters contribute at most

    s²/[π²(1/4+s²)(1/256+s²)].

Positive two-column coarea and Q4 scaling yield

    P_s≤R2,min(k,m) s^−2/[16π⁶(1/4+s²)(1/256+s²)]
       ≤Dbig,k s^−2,   Dbig,k=64R2,min(k,m)/π⁶.           (4)

The root J2 cancels only because the frozen R2 is the normalized Q4/J2 coefficient. It is not silently set to one. The three exact R2 entries and their Gram/kappa definitions are read from033; no new peak computation occurs.

## 4. Non-tight deletion, missing parameter and complete bad mass

The original actual parameter events are the disjoint exhaustive sets

    E: |a_v|≥1/2 and |a_w|≥1/16,
    E_v: |a_v|<1/2,
    E_w: |a_v|≥1/2 and |a_w|<1/16.

Let P_s be the density for E and B_s the combined actual complement. They sum to p_s and separately satisfy0≤P_s,B_s≤p_s. For the pointwise upper estimate only, E_w may be enlarged by dropping its v restriction. The enlarged events need not be disjoint and their sum is not claimed≤p_s.

Here is the positive deletion argument without a fictitious tight remaining frame. Fix the retained small b=v or w. After r=0,1,2 non-b deletions and orthogonal projection P, the remaining31−r columns have Gram

    32I−(Pb)(Pb)^T,

between26I and32I because |Pb|²≤6; already deleted projected columns vanish. For any such matrix C and parameters a, |Ca|²≤32Σa_i². For their product F_s,

    (|Ca|²+32Ns²)F_s ≤32(s/π)Σ_i ∏_(j≠i)f_s(a_j),       (5)

where N is the actual number of remaining columns. Positive pushforward and line coarea give the deletion inequality. Consequently the three successive denominators have precisely32·31s²,32·30s²,32·29s², namely992,960,928. No tight31I approximation or artificial parameter-time change is used. Products of successive projected-column norms give the true J3. Three distinct body lines are independent, as established in043, so those factors are nonzero. The Gram lower bound also supplies continuous quotient densities.

After the three deletions, the dimension3 quotient still omits b. On this quotient its symbol is ψ−|b·ξ|. The frozen043 comparison ψ≥10||ξ||1 and |b·ξ|≤||ξ||1 give ψ−|b·ξ|≥(9/10)ψ. Fourier peak comparison and three-dimensional homogeneous change of variables therefore require the factor(10/9)^3. This only compares peaks; it does not change the common time in the original event law.

On the actual small-a_b domain, cap the first two denominators using(2). Only then retain the third integral and extend this positive integral to all real a_b. The043 exact Cauchy average gives

    ∫_R f_s(a)/(|z−a b'|²+928s²) da
       ≤(1+|b'|/sqrt(928))/|z|².

Use original z=P_(span(i,j))⊥x with |z|≥g_k from(3), and |b'|≤sqrt(6). Isolated shifted zero distances are included by this formula; none is discarded. The first-two-gap estimates are not being asserted outside the original small-parameter domain during this extension.

The complete frozen ordered-triple constant for one omitted body line is C_rem=44743079/(2970π⁶), including all26970 ordered triples and their true J3 factors. Signed-coordinate symmetry makes it the same for either retained body line, independently of k. Thus

    L_k=(10/9)^3 C_rem(1+sqrt(6/928))/g_k²,
    Mv,k=L_k/(hv0²hv1²), Mw,k=L_k/(hw0²hw1²),
    B_s≤M_k=Mv,k+Mw,k.                                   (6)

This keeps every actual complementary event, all other large parameters and all cancellations already paid by043. The remaining30-column frame for the two-large branch has Gram32I−vvᵀ−wwᵀ≥20I, so that branch's continuous pointwise argument remains valid as well.

## 5. Joint marginal for every B size, including m=1

Choose one fixed j∈B. The exact body32 marginal and remaining-four-coordinate positive orthant factor are those independently reviewed in049: (1/2)f_(16s)((1+z)/2)f_(16s)((1−z)/2), and1/16 respectively. Symmetry is applied only to full p_s.

On the zero-sum B block, the coordinate functional has norm sqrt((m−1)/m), equal to zero for m=1. Since the actual residual has norm<rho, its retained coordinate lies in

    L_m=1/4−sqrt((m−1)/m)rho,
    U_m=3/4+sqrt((m−1)/m)rho.

Both lie strictly between0 and1. The singleton case has exactly[1/4,3/4], not an empty quotient or a spurious residual factor. Positivity and extension of the other four face coordinates to their positive orthant give

    ∫_DA p_s dy≤C_m s²,
    C_m=(128/π²)∫_(L_m)^(U_m)(1−z²)^−2dz.                 (7)

There is no multiplicity m because j is fixed. Apply domination P_s,B_s≤p_s to(7), not branch reflection symmetry. Interpolation at q=6/5 and positive time Minkowski/monotone convergence then give

    ||S_epsilon||_(Lq(DA))≤U_k,
    U_k=(3/7)Dbig,k^(1/6)C_m^(5/6)epsilon^(7/3)
        +(3/8)M_k^(1/6)C_m^(5/6)epsilon^(8/3).             (8)

The k1 and k5 cases require no alteration of the original probability laws or events.

## 6. k3 source reuse and overlap-safe31-domain union

For k=m=3 the geometry and formulas are the already accepted043/049 ones. This review does not evaluate those scalar formulas again. It reads the full exact rational total upper endpoint from its own frozen049 independent result; the author's k3 total-budget entry is additionally checked to be the saved049 total-budget entry.

For each k the permutations of the five free coordinates are transitive on the A subsets of that size and preserve the model, domain and face measure. Hence n_k=binom(5,k−1)=(1,5,10,10,5). Their sum is31. These weights cannot be merged by k↔6−k, because coordinate1 remains in the distinguished A block.

For the same nonnegative function S_epsilon and q=6/5,

    ∫_(union_A DA) S_epsilon^q
      ≤Σ_A∫_DA S_epsilon^q≤Σ_k n_k U_k^q.

Taking the1/q power gives the claimed union bound. The inequality permits every overlap and never calls the domains disjoint. This is only the union of these specified middle domains.

## 7. Independent numerical acceptance and canonical floor correction

Exactly six new relation/scalar configurations were run: one common symbolic structure group; four Arb768 cases k1,2,4,5; one weighted union. Gaps were factored as sqrt(k)[sqrt(m/6)/4−rho] or /8−rho, the primitive evaluated using atanh, and new short terms as sixth roots of Dbig C_m^5 epsilon^14 or M_k C_m^5 epsilon^16. Each new coefficient/term interval is compatible with the author's complete saved dyadic certificate. The union uses exact upper endpoints of our four new cases and the frozen049 k3 endpoint, so all input directions are upper-bound directions.

Its independent outward upper endpoint satisfies strictly

    union_upper <4689838525803/10^18.                     (9)

The exact positive margin and every new upper endpoint are retained in results.json. No printed decimal determines acceptance.

The author's first rational-only saved-data checker failed at a strict floor sanity assertion before its case checks. Its canonical correction to non-strict≤ is valid: floor(T·π_lower)/T is a lower bound for π and may exactly equal the chosen rational floor constant. Equality does not reverse any π lower-bound direction. Using a positive lower π in denominators gives upper bounds, as required. Its positive atanh-series remainder also has the correct upper direction, while primitive(U)_upper−primitive(L)_lower enlarges the integral. This review did not execute either version of that checker; independent Arb endpoint acceptance in(9) stands on its own. The failed first receipt is retained by the author and not hidden.

No core gap was found. All tails above epsilon, spatial mean endpoints outside[1/4,3/4], other regions of the face, full I, global g and Cn remain unpaid.
