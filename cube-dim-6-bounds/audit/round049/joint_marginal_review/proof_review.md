# A049 joint marginal: independent proof review

This review fixes the original32 body columns (1,s2,…,s6), original common semigroup time s, A={1,2,3}, B={4,5,6}, U0=span(1_A,1_B), rho=1/64 and epsilon=1/1000. The face measure is dy in x=(1,y). The exact domain is

    D={y∈[0,1]^5 : dist((1,y),U0)<rho,
                         mean_B((1,y))∈[1/4,3/4]}.

The F043-reviewed pointwise constants and event account are frozen inputs, read from their complete saved intervals. No043 pointwise formula or old short-budget formula is recomputed.

## 1. Exact retained marginal and its normalization

Fix j=4. Among the body columns, exactly16 have j-sign +1 and16 have j-sign −1, since four signs remain free. Let U and V be sums of the independent standard Cauchy(s) variables in these two disjoint sets. Their characteristic functions are exp(−16s|η|), so U,V are independent Cauchy(16s). The retained random vector is

    (X1,X4)=(U+V,U−V).

The forward determinant has absolute value2, so the inverse density Jacobian is1/2. At (1,z), direct substitution gives

    p2_s(1,z)
      =(1/2)f_(16s)((1+z)/2) f_(16s)((1−z)/2)
      =2048s²/[π²(((1+z)²+1024s²)((1−z)²+1024s²))].          (1)

There is no additional factor for an unoriented frame: the32 columns already have their stated independent parameters.

The slice integration identity is valid at these fixed retained coordinates, not only at almost every pair. Choose basis columns ones, flip4, flip2, flip3, flip5, flip6. Their determinant has absolute value32. Order the physical rows as1,4,2,3,5,6. The retained2×2 pivot determinant has absolute value2; the four-coordinate Schur map has absolute determinant16. In the ordinary positive full-rank basis integral for p_s, integration over the four other physical coordinates changes those four free basis variables with Jacobian16. It cancels the original1/32 to leave1/2. Tonelli over the remaining30 variables then yields the two16-fold Cauchy convolutions in(1), pointwise in the two retained coordinates. This uses positive integrals and an invertible change of variables, not conditioning on a null fiber.

## 2. The four-coordinate orthant factor belongs to full p only

Flip any coordinate in {2,3,5,6}. The retained coordinates1,4 stay fixed, and sign-bit reversal permutes the body column set. Independence of identical Cauchy variables therefore preserves the full law, hence the full continuous density p_s. The16 orthants of the four integrated coordinates have equal integrals. Their boundaries are Lebesgue-null in four dimensions, so

    ∫_(remaining four≥0) p_s(1,z,rest)drest=p2_s(1,z)/16.       (2)

Neither retained coordinate contributes a further factor. In particular1/32 would be wrong here. A distinguished-parameter event is generally changed by these flips, so no such equality is assigned to an event-restricted density.

## 3. Domain inclusion and the joint L1 bound

Write x_B=u1_B+q_B, with u its actual mean and q_B in E_B={q∈R³:Σq=0}. This residual space has dimension2. Orthogonal projection gives dist(x,U0)²=|e_A|²+|q_B|², so |q_B|<rho on D.

The coordinate functional q↦q4 on E_B is represented by e4−1_B/3, of squared norm2/3. Thus for the same fixed j=4, at every point of D,

    L≤x4≤U,
    L=1/4−sqrt(2/3)/64, U=3/4+sqrt(2/3)/64.                 (3)

The assertion holds for each coordinate in B, so a fixed j is sufficient; there is no adaptive-coordinate union or factor3. Since sqrt(2/3)<1, one has0<L<U<1.

All four remaining face coordinates are nonnegative. Extend their upper box limits from1 to infinity, keeping x4 in[L,U]. Positivity, (1) and(2) imply

    ∫_D p_s dy
      ≤(1/16)∫_L^U p2_s(1,z)dz
      ≤ Cjoint s²,
    Cjoint=(128/π²)∫_L^U(1−z²)^−2 dz.                       (4)

The last inequality drops positive1024s² terms in the denominators; (1−z²)² is positive throughout this exact interval. All measures are the original face coordinate measures. No residual Jacobian is inserted into (4), because residual coordinates were used only to prove set inclusion.

An independent primitive is

    F(z)=[z/(1−z²)+atanh(z)]/2,
    F'(z)=(1−z²)^−2,   0<z<1.

This is equivalent to the author's log-ratio primitive. The review evaluates the atanh form and obtains an outward Arb768 enclosure for the exact Cjoint.

## 4. Complete actual-event bookkeeping

With v=ones and w=(+++---), the actual three sets are

    E: |a_v|≥1/2 and |a_w|≥1/16,
    E_v: |a_v|<1/2,
    E_w: |a_v|≥1/2 and |a_w|<1/16.

They are disjoint and exhaustive, including their chosen threshold boundaries. All other parameters are unrestricted. If P_s is the pushforward density of E and B_s the sum for E_v,E_w, then P_s+B_s=p_s and separately0≤P_s,B_s≤p_s. This is the actual complement, not the sum of two enlarged events.

The F043 proof bounds these same densities on this same D by

    P_s≤Dbig s^−2,   B_s≤Mbad.

Here Dbig is043 configurations[3].C2 and Mbad is043 configurations[4].M_combined. Its proof may enlarge E_w to {|a_w|<1/16} while bounding the supremum; that enlarged event may overlap E_v. Such overlap does not invalidate the pointwise sum upper bound, but it would invalidate a claim that the enlarged sum is≤p_s. The new interpolation instead uses actual B_s≤p_s, hence by(4)

    ||P_s||_1(D)≤Cjoint s²,   ||B_s||_1(D)≤Cjoint s².          (5)

No branch reflection symmetry is used. Actual B receives one combined L1 bound, not two separate whole-density L1 charges. Every parameter event already covered by043 remains included.

## 5. Interpolation and positive-time integration

For nonnegative f, ||f||_(6/5)≤||f||_∞^(1/6)||f||_1^(5/6). Combining the pointwise constants with(5) yields

    ||P_s||_(6/5)(D)≤Dbig^(1/6)Cjoint^(5/6)s^(4/3),
    ||B_s||_(6/5)(D)≤Mbad^(1/6)Cjoint^(5/6)s^(5/3).

For positive lower time cutoffs apply Minkowski and the triangle inequality. Then let the cutoff decrease to zero using monotone convergence of the nonnegative time integral and its6/5 power. The result is

    ||∫_0^epsilon p_s ds||_(6/5)(D)
      ≤(3/7)Dbig^(1/6)Cjoint^(5/6)epsilon^(7/3)
        +(3/8)Mbad^(1/6)Cjoint^(5/6)epsilon^(8/3).             (6)

The powers are exactly −2/6+2·5/6=4/3 and2·5/6=5/3; integrating adds1 and gives prefactors3/7 and3/8. No regularity of full G at singular strata is assumed. Singular infinite point values on null sets do not break this Lq argument.

## 6. Independent interval acceptance

After final parent pins matched, the full author whitelist and frozen043/F043 inputs were hash-checked. Only the complete binary mantissas/exponents of Dbig,Mbad and the old budget are imported. They become exact Fractions before any new arithmetic; no printed decimal is reconstructed.

The independent new terms are evaluated as

    (3/7) sixthroot(Dbig Cjoint⁵ epsilon¹⁴),
    (3/8) sixthroot(Mbad Cjoint⁵ epsilon¹⁶).

For each, the saved constant lower/upper endpoints are evaluated separately using positivity and monotonicity; the sixth root is computed as exp(log(value)/6) with outward Arb768. The total upper is the exact rational sum of the two independent upper endpoints. It satisfies

    independent_total_upper <66002272277/250000000000000000.  (7)

The strict positive margin is retained as a Fraction in results.json. Author interval overlap is additionally checked, but acceptance rests on our own strict upper endpoint comparison. The saved old-budget lower endpoint is above the new rational threshold, a strict improvement of upper-budget formulas. It is not a lower bound on the actual norm and is not an actual norm reduction factor.

No old043 formula was reevaluated, no author writer/checker was executed, and no kernel or full enumeration was used. This review certifies precisely the same single middle D and0<s≤.001. Long time, other planes, endpoints, whole-face I, global g and C6 remain unpaid.
