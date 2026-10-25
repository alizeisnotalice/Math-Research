# Full-time S1 cap from the two actual parameter events

The frozen A041 positive decomposition and F041 PASS have been pinned in dependency_release.json before the sole scalar batch. It does not introduce a time-averaged model or a new spatial point. Fix the original body32 Cauchy model, the positive face F={(1,y):y∈[0,1]^5}, q=6/5 and rho=1/64. Let C={x∈F:dist(x,S1)<rho}. By034, the all-ones line attains the nearest S1 distance on F, and

    V=|C|=(4√6 π²/45)rho^5.

This is the five-dimensional face-coordinate Lebesgue volume, with no extra31 or32 multiplier.

## 1. The two events and the full time partition

For v0=1, A041 splits the actual Cauchy parameter a_v0 at |a_v0|=1/2:

    p_t=P_t+B_t,  P_t≥0, B_t≥0,

where P is the |a_v0|≥1/2 event and B its complement. The other31 parameters are unrestricted in both events. Multiple large jumps, other large directions and cancellation events are included, not omitted. Assume the A041 all-t bounds on C:

    P_t ≤ D t^-4/(1/4+t²),  D=R1/(32π^6),
    B_t ≤ Crem/∏_(i=0)^2(c_i²+A_i t²),

    R1=67320683/3251298631680, Crem=44743079/(2970π^6),
    (A0,A1,A2)=(992,960,928),
    c0=√6/2−√6rho,
    c1=√(10/3)/2−√6rho,
    c2=1/√2−√6rho.

All three c_i are positive at the fixed rho. Put Cbig=4D. A041/B036 also gives by P_t≤p_t

    ∫_C P_t dy ≤ t/[π(1+1024t²)] ≤ t/π.

No orthant symmetry is asserted for P_t itself. Interpolating this L1 domination with P_t≤Cbig t^-4 yields

    ||P_t||_q ≤ Cbig^(1/6) π^-5/6 t^1/6.

The whole-density slice bound can also be checked directly: the first coordinate is the sum of32 independent Cauchy(t) variables, so its marginal density at1 is32t/[π(1+1024t²)]. Coordinate sign flips preserve the full body32 law, making the five other coordinates even. The positive orthant therefore contributes exactly1/32 of that marginal; C is a subset of it and P_t≤p_t. No analogous symmetry of the selected event is required.

For any epsilon>0,

    ||∫_0^epsilon P_t dt||_q
       ≤S(epsilon):=(6/7)Cbig^(1/6)π^-5/6 epsilon^7/6.     (1)

For the entire remaining P tail, use its pointwise bound and the exact cap volume:

    ||∫_epsilon^∞ P_t dt||_q
       ≤V^(5/6) D R(epsilon),                            (2)

    R(epsilon)=∫_epsilon^∞ t^-4/(1/4+t²)dt
       =4/(3epsilon³)−16/epsilon+32 atan[1/(2epsilon)].    (3)

To prove (3), the integrand is exactly

    4t^-4−16t^-2+16/(t²+1/4).

Differentiate the right side of(3) to obtain minus the integrand and note its zero limit at infinity. The arctangent branch is positive since epsilon>0; equivalently atan(1/(2epsilon))=π/2−atan(2epsilon).

We choose only one epsilon, using R≤4/(3epsilon³). The loose bound is

    S(epsilon)+V^(5/6) Cbig/(3epsilon³).

Its derivative vanishes precisely when epsilon^(25/6)=(Cbig πV)^(5/6), so the unique positive minimum occurs at

    epsilon*=(Cbig πV)^(1/5).                            (4)

The derivative changes from negative to positive there, and the loose bound tends to infinity at either endpoint. We evaluate (1)-(3) only at this single analytic value; we do not optimize the exact tail, scan epsilon or alter the Cauchy model. In particular epsilon* is an auxiliary split of the already fixed Green time integral, not a new time-average parameter.

## 2. The complementary branch on the whole positive time axis

For positive a,b,c, set J(a,b,c)=∫_0^∞∏(t²+a_i²)^-1 dt. First, the elementary two-factor partial fraction identity gives

    I2(a,c)=∫_0^∞[(t²+a²)(t²+c²)]^-1 dt
            =π/[2ac(a+c)].

For a≠b, J=[I2(a,c)−I2(b,c)]/(b²−a²). Exact cancellation gives the positive symmetric formula

    J(a,b,c)=π(a+b+c)/[2abc(a+b)(a+c)(b+c)].              (5)

This derivation fixes both the sign and every Fourier-free factor. On any compact subset of positive parameter triples the original integrand is dominated by(t²+m²)^-3 for some m>0. Dominated convergence extends(5) to every repeated-root case. Independently,

    J(a,a,c)=−(1/(2a))∂_a I2(a,c)
             =π(2a+c)/[4a³c(a+c)²],
    J(a,a,a)=3π/(16a^5).

Thus no divided difference with a small or zero root separation is needed in the final evaluation. Zero a_i is not admitted; positivity of c_i is essential.

Set a_i=c_i/√A_i. Pulling out ∏A_i gives

    JB=∫_0^∞ [∏(c_i²+A_i t²)]^-1 dt
       =J(a0,a1,a2)/(A0 A1 A2).                         (6)

Since the bound is uniform on C, the **entire** complementary branch is paid by

    ||∫_0^∞ B_t dt||_q ≤ V^(5/6) Crem JB.                (7)

## 3. Complete cap budget and overlap ledger

For finite positive time truncations, Minkowski applies to the three nonnegative pieces P on(0,epsilon*), P on[epsilon*,∞), and B on(0,∞). Then monotone convergence extends to the full Green function and proves

    ||G||_Lq(C) ≤ S(epsilon*)
                  +V^(5/6)D R(epsilon*)
                  +V^(5/6)Crem JB =:Ucap.               (8)

The two parameter events are disjoint and exhaustive. The two P time ranges are disjoint and exhaustive; B is included on both time ranges through its one whole-axis integral. We do not assume the resulting spatial functions have disjoint supports. Minkowski pays their sum, so there is no negative overlap correction or lost event. The Cauchy threshold boundary and the spatial cap boundary have measure zero and do not change the integral.

This one cap contains its S1 corner and every incident S2 crossing; these are already covered by the actual p=P+B inequality and do not require another plane multiplier. When later joining a tube budget, either use a positive overlapping union bound or explicitly partition its complement by subtracting C. This theorem does not treat any region outside C. In particular

    ∫_C G^q ≤ Ucap^q

is only a cap contribution, not the complete face I. Neither the separate original R_c floor g≥1/65000 nor Hnorm≤1/14000 has been proved here. C6 retains its established upper bound64.
