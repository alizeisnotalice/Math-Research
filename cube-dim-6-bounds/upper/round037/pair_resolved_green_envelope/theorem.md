# Pair-resolved positive deletion and an integrable Green envelope

Use the original body32 Cauchy frame V={(1,±1,...,±1)} and orthogonal quotient Lebesgue measures. All distances below are actual Euclidean distances to the indicated individual subspace. The B008 positive deletion lemma and the exact033 peaks, now independently accepted by F035, are the inputs. Write

    r=|x|, d_v=dist(x,span v), ρ_vw=dist(x,span(v,w)),
    D0=r²+1024t², Dv=d_v²+992t², Dvw=ρ_vw²+960t².

## 1. Keep the positive branches before taking any common distance

B008 proves the genuine quotient inequality

    q_U,t(P_Uperp x)
      ≤32t/[π(dist(x,U)²+32(32−dim U)t²)]
        Σ_(i not deleted) |P_Uperp v_i|^-1 q_(U+span v_i),t.

For the empty U the distance is r. For one or two deleted directions all remaining projected vectors are nonzero. Every three distinct body lines are independent. The product of sequential projection lengths is J_k=√det Gram of the deleted directions. This is the coarea factor; no common coordinate chart or additional Jacobian is introduced here.

Let Q4(v,w) and Q3(v,w,u) be exact quotient peaks at time1. After two deletions, bound the actual q_(v,w),t by BOTH its peak Q4 t^-4 and its positive third-deletion estimate. Taking their minimum is valid separately in each positive summand. This gives

    p_t(x) ≤ Σ_(ordered v≠w) [D0 Dv]^-1
                min{A_vw t^-2, C_vw/Dvw},                    (1)

where

    A_vw=(32/π)² Q4(v,w)/J2(v,w),
    C_vw=(32/π)³ Σ_(u≠v,w) Q3(v,w,u)/J3(v,w,u).             (2)

Every ordered pair and every third direction occurs. This is an upper estimate of a positive sum, not a cancellation or an allocation of each sample to one nearest plane.

033 gives Q_(6−k)/J_k=R_k/(2π)^(6−k), so the Gram square roots cancel **only after the correct orthogonal coarea factor has been retained**. Thus, for pair orbit k=min(Hamming(v,w),6−Hamming(v,w)),

    A_k=64 R2,k/π^6,     C_k=4096 S_k/π^6,
    S_k=Σ_(u≠v,w) R3,class(v,w,u).                          (3)

No new frame enumeration is needed to obtain these small sums. Normalize the first two rows to agreement/difference blocks of sizes k and6−k. Choose the third row's signs within each block, and divide all64 signed choices by global negation; discard exactly the four dependent signed choices ±v,±w. Binomial block counts give

    S1=10 R4110+20 R3210,
    S2=2 R4110+8 R3210+6 R2220+8 R3111+6 R2211,
    S3=12 R3210+18 R2211.                                 (4)

Each coefficient list sums to30, the exact number of third body directions. For example blocks3+3 are either constant (2 sign choices) or split2+1 (6 choices). One constant and one split block gives24 signed choices, hence12 of class3210; two split blocks give36, hence18 of class2211. For blocks2+4 the nondependent signed counts are respectively4,16,12,16,12 for classes4110,3210,2220,3111,2211, before dividing by2. This accounts for all mixed sign rows.

The actual ordered-pair orbit counts are192,480,320: for each of32 first directions the partner has binomial(5,h) choices, h=1,...,5, and h is folded with6−h. These counts sum to992. They are counts in (1), not unordered counts.

The Q5 input is also retained: first deletion alone gives

    p_t(x) ≤ [32 R1/π^6] t^-4/D0.                         (5)

One may take the minimum of the complete right sides of (1) and (5). Dropping that optional improvement below leaves a valid, simpler envelope.

## 2. Closed Green envelope and its integrability

For A,C>0 set

    M(A,C)=∫_0^∞ min{A/s²,C/(1+960s²)} ds.

If C≤960A, M=Cπ/(2√960). If C>960A, let s0=√[A/(C−960A)]; then

    M=C/√960 atan(√960 s0)+A/s0.                          (6)

All three actual pair types satisfy C>960A, checked as the exact rational inequality 4096S_k>960·64R2,k, before Arb is used. Substituting t=ρ_vw s and using D0≥r², Dv≥d_v² in (1) proves, for x outside S2,

    G(x) ≤ Φ(x):=Σ_(ordered v≠w)
             M_k/[r² d_v² ρ_vw].                         (7)

At zero denominators define Φ=+∞, preserving the upper bound. Formula (7) is a global **upper envelope**, not an asymptotic equality for the true G. In particular its coefficients are not asserted to be the actual two-jump coefficients of G.

On the positive face F={(1,y):0≤y_i≤1}, r≥1. Each individual summand belongs to L^q(F), q=6/5. To see this, only the all-ones body line meets F. If v is any other body representative, compactness gives d_v bounded below on F, and the transverse distance to the affine pair intersection has codimension4 in the face; ρ^-q is locally integrable for q<4. More explicitly, the map P_(span(v,w))perp on the five-dimensional tangent hyperplane has rank4: its one-dimensional kernel is span(v,w)∩{h1=0}.

For v=1, put x=v+h, h1=0. Then d_v≥|h|/√6. Decompose h=sℓ+z, where ℓ spans span(v,w)∩{h1=0} and z is its four-dimensional orthogonal complement inside the face tangent. The nonzero singular values of the projection to span(v,w)^perp are at least1/√6: the smallest squared value is |P_span(v,w)e1|²≥1/6 since v belongs to that plane. Hence ρ_vw≥|z|/√6. The q-th power is dominated by a constant times

    |h|^-2q |z|^-q.

In five-dimensional polar coordinates its radial integral has exponent4−3q and its angular factor near z=0 has exponent3−q. They are integrable when3q<5 and q<4, both true for6/5. A finite sum preserves L^q. Intersections between different planes are therefore covered without assuming disjointness or omitting mixed terms. This proves existence of a legal finite L^q majorant; it does not evaluate its whole-face norm.

## 3. One unsigned plane: only two ordered terms are singular in its middle

Fix A={1,2,3}, B={4,5,6}, v=1, w=1_A−1_B, U=span(v,w), and

    x(u)=1_A+u1_B, 1/4≤u≤3/4.

The individual distances are exactly

    r²=3(1+u²), d_v²=(3/2)(1−u)², d_w²=(3/2)(1+u)².        (8)

A pair plane contains exactly its two body lines: if av+bw is another cube vertex, the nonempty agreement and difference blocks force |a+b|=|a−b|=1, hence ab=0. Thus each plane appears exactly twice in the ordered pair sum; it cannot be silently counted once. At this positive two-level middle point, U is the only pair plane containing x(u).

Here is a uniform quantitative argument that also retains every other plane. Any other unsigned partition has a group containing coordinates with the two distinct values1 andu; its distance is at least(1−u)/√2. For a signed plane with the same partition A/B but a genuinely different span, an opposite sign pair in a block gives distance at least√2 for the A block, or √2u for B. Other signed partitions dominate their corresponding unsigned distance by034. Consequently **all other495 planes** have distance at least

    min{(1−u)/√2,√2u} ≥1/(4√2).                          (9)

The all-ones line is the nearest body line on the positive face, and its middle distance is at least√(3/32). The corresponding statements at u=0 or u=1 do not hold: mixed signed planes can meet the u=0 endpoint, while d_v tends to0 at the S1 endpoint u=1. Those endpoints require their own analysis; they are not paid by the middle estimate.

For approaches to x(u) off U, with ρ=dist(x,U), (9) makes all990 other ordered terms bounded. Therefore the upper envelope itself has the exact local limit

    lim_(x→x(u)) ρ Φ(x)
      =LΦ(u)=M3/r²(1/d_v²+1/d_w²)
      =4M3/[9(1−u²)²].                                  (10)

This is the **actual coefficient of Φ**, with true ordered distances and both orders included. It is NOT a claim that ρG tends to this number. The common-distance collapse of the same M constants would instead have coefficient

    Lcollapsed(u)=Mtotal/[r² d_v²],
    Mtotal=192M1+480M2+320M3.

At the three fixed u values1/4,1/2,3/4 the ratio Lcollapsed/LΦ is about791.630,968.955,1055.084. This large change comes from localizing the singularity to its two genuine terms, leaving990 bounded terms, and retaining the large d_w distance in the reversed order. It is additional information beyond replacing the old global B3 by better peaks.

## 4. One finite middle tube, with all the background paid

Use only the fixed radius δ0=1/64. Let D consist of positive face points with dist(x,U)=ρ<δ0 and mean_B x=u∈[1/4,3/4]. In the034 coordinates x_A=α1_A+e, x_B=u1_B+z, with Σe=Σz=0 and e1=1−α,

    |x−x(u)|²=ρ²+3e1²≤3ρ²,

since |e1|²≤(2/3)|e|². The1-Lipschitz property of subspace distance yields uniform positive gaps

    a=√(3/32)−√3δ0>0  for every body line,
    g=1/(4√2)−√3δ0>0 for every pair plane other than U.

Thus (7), with all992 terms, gives

    G≤Φ≤Atube/ρ+Bbackground,
    Atube=2M3/a²,
    Bbackground=(Mtotal−2M3)/(a²g).                       (11)

No other-plane term has been dropped; mixed planes and overlapping tubes are included in the990-term background. There is no probability or independence assumption about planes.

034 gives face Jacobian √(3·3)=3 and an exact1/3 cone fraction for radial integration, with middle length1/2. Extending the remaining box constraints positively therefore gives

    |D|≤(v4/2)δ0^4,
    ∫_D ρ^-q ≤[2v4/(4−q)]δ0^(4−q), v4=π²/2.

Minkowski now proves the fully paid single-tube estimate

    ||G||_Lq(D)
      ≤Atube [2v4/(4−q)]^(1/q) δ0^(4/q−1)
       +Bbackground [(v4/2)δ0^4]^(1/q).                  (12)

The sole registered tube diagnostic gives (12) about8.00025266411e-5, approximately0.9000284 times the entire old target1/11250. This is only one middle tube, not its endpoints, other plane types, their union or the whole face. It is a useful but currently expensive positive budget, not a new I or C6 certificate. No extrapolation or sum over uncomputed tubes is claimed.
