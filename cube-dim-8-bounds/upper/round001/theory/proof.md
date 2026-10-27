# R001 — Centered cubes: a self-contained cap upper bound in every dimension

Date: 2026-10-27. Author role: Sol theory reconstruction. Status: the analytic argument below is complete as written; fresh project review remains the parent's responsibility. No old PASS verdict is used as a premise. No numerical upper value, optimized parameter, published novelty claim, or variable-background extension is asserted.

## 1. Operator, theorem, and scope

Write Q(y,r)=y+(-r,r)^n, where n is a positive integer and r>0 is **one common radius for all coordinates**. All volumes and all integrals in this note use Lebesgue measure on R^n. Define

    M_n f(x)=sup_{r>0} (2r)^(-n) ∫_{Q(x,r)} |f(y)| dy,
    C_n=sup_{f∈L¹, ||f||₁>0, λ>0} λ|{M_n f>λ}|/||f||₁.

Thus C_n is the full weak-(1,1) constant, using arbitrary real or complex L¹ inputs, all centers, and every positive radius. We will prove that, for every t>1,

    C_n ≤ [(1+1/t)^(n/2)+sqrt(B_n(t))]²,                 (1.1)

where

    B_n(t)=1-[(1+1/t)/2]^n+n log t
           -(n/2) ∫₁ᵗ [1+1/(2t)-1/(2ρ)]^(n-1) dρ/ρ.   (1.2)

The square root is well defined because B_n(t)>0, proved below. In particular,

    C_8 ≤ [(1+1/t)^4+sqrt(B_8(t))]²,
    B_8(t)=1-[(1+1/t)/2]^8+8 log t
           -4 ∫₁ᵗ [1+1/(2t)-1/(2ρ)]^7 dρ/ρ.          (1.3)

No restriction of the maximal operator to a finite test configuration occurs in (1.1). Finite families appear only inside a covering argument applied to arbitrary compact subsets of a true operator level set.

The ambient background is the constant Lebesgue density 1. This is used in the cube volumes, in the masses of the source increments, in dilation and translation changes of variables, and in the shell calculation. The input |f| itself need not be constant, positive everywhere, bounded, compactly supported, a product, or a finite source. The geometric proof is not a bound for maximal operators normalized by arbitrary variable measures.

## 2. Keep a finite enlargement parameter until the operator bound is proved

Fix t>1, a>2t, and τ>0. Set

    q=1+1/t,             b=a q,
    α=a/(a+1),           c=(b-1)/(a+1)>1,
    g_t(ρ)=1+1/(2t)-1/(2ρ),
    I_n(t)=∫₁ᵗ g_t(ρ)^(n-1) dρ/ρ,
    H_n(a,t)=1-α^n(q/2)^n+n log t
             -(n/2)α^(n-1) I_n(t),
    R_n(a,t)=a^(-n)+α^(-n)H_n(a,t),
    A_n(a,t)=(q+1/a)^n.

We first prove the genuine full-operator inequality

    C_n ≤ [1+A_n(a,t)/τ][τ+R_n(a,t)].                    (2.1)

The self term a^(-n) in R_n is essential at every finite a. No proof step sets a=∞ or drops that term. Only after (2.1) is established for all admissible finite a will its scalar right side be sent to its limit.

## 3. Selection from an arbitrary finite family

Start with any finite family of open cubes. Repeatedly select an eligible cube of maximal radius among all unselected cubes. Once Q_j=Q(y_j,r_j) has been selected, freeze

    D_j=Q_j \ ⋃_{i<j}Q_i,
    v_j=|D_j|,
    w_j=v_j/(2a r_j)^n.

At any later stage, a still unselected cube centered at y is eligible precisely when

    ∑_{j already selected} w_j 1_{bQ_j}(y) ≤ τ.           (3.1)

Here dQ means the concentric cube with radius d times the original radius. Stop when no eligible cube remains. There are only finitely many cubes, so this stops. Zero increments and equal radii are allowed. All weights are nonnegative; once frozen, they never change.

The D_j are pairwise disjoint measurable sets, and

    ∑_j v_j=|⋃_{j selected} Q_j|.                        (3.2)

Since eligibility can only be lost, the selected radii are nonincreasing. Indeed any cube eligible at a later stage was eligible at the preceding selection; maximality of the preceding chosen radius then gives the ordering.

For each unselected cube Q_β=Q(y_β,r_β),

    ∑_{j:r_j≥r_β} w_j 1_{bQ_j}(y_β)>τ.                  (3.3)

For completeness, if some selected radius is smaller than r_β, consider the first such selection. At that time Q_β was unselected and could not have been eligible, because it has the larger radius. All earlier selected radii are at least r_β. Their rejection sum was already >τ. If no selected radius is smaller than r_β, at termination Q_β is ineligible and every selected radius is at least r_β. Either case proves (3.3).

For x∈Q_β and every contributing j in (3.3),

    ||x-y_j||∞ ≤ ||x-y_β||∞+||y_β-y_j||∞
               < r_β+b r_j ≤ (b+1)r_j.

Consequently the final function F(x)=∑_j w_j1_{(b+1)Q_j}(x) is >τ at every point of the union of unselected cubes. Integrate this inequality, then add the selected union and use (3.2):

    |⋃_{all supplied}Q_β|
      ≤ ∑_jv_j + τ^(-1)∫F(x)dx
      = [1+((b+1)/a)^n/τ]∑_jv_j
      = [1+A_n(a,t)/τ]∑_jv_j.                            (3.4)

The dilation is b+1, including the radius of the rejected cube. No star-cube disjointness was assumed.

## 4. Pointwise overlap: earlier terms, ties, and self payment

For this selected family let

    W(x)=∑_j w_j1_{aQ_j}(x).

We claim W(x)≤τ+R_n(a,t) for every x. If W(x)=0 it is immediate. Otherwise, among all **positive contributing terms** choose k of minimum radius, and among terms at that minimum radius choose the **last selected index**. Translate x to 0 and divide all coordinates and radii by r_k. Incidences and weights are unchanged, so now r_k=1 and 0∈aQ_k.

Every other positive contributor has r_j≥1 and index j<k. To see the latter, if j>k, nonincreasing selected radii give r_j≤r_k=1. Hence equality holds, contradicting the last-index choice. This is the precise treatment of equal-radius ties. Terms with w_j=0 require no payment.

The self term, absent from the earlier eligibility test at y_k, satisfies separately

    w_k ≤ |Q_k|/|aQ_k|=a^(-n).                           (4.1)

Call another positive contributor captured if y_k∈bQ_j. All such contributors have j<k, and they form a subset of the nonnegative sum in the actual test (3.1) when Q_k was selected. Their combined weight is therefore ≤τ.

An uncaptured positive contributor satisfies

    ||y_k-y_j||∞≥b r_j,
    ||y_j||∞<a r_j,
    ||y_k||∞<a.

If r_j≥t, the latter two inequalities give

    ||y_k-y_j||∞<a+a r_j≤a(1+1/t)r_j=b r_j,

a contradiction. Thus every unpaid positive term has

    1≤r_j<t.                                             (4.2)

This includes earlier terms with radius exactly 1; all of them really are earlier because of the tie rule.

## 5. Turn unpaid increment masses into one geometric cap integral

For y∈D_j of an uncaptured positive contributor, use D_j⊂Q_j to get

    ||y_k-y||∞ ≥ ||y_k-y_j||∞-||y_j-y||∞
                >(b-1)r_j.                              (5.1)

Set z=y/(a+1), u=y_k/(a+1), and m(z)=max(1,||z||∞). Since 0∈aQ_j,

    ||u||∞<α,
    ||z||∞<r_j<t,
    r_j≥m(z),
    ||u-z||∞>c r_j≥c m(z).                              (5.2)

Write E_j=(a+1)^(-1)D_j. These E_j remain pairwise disjoint. The larger common comparison domain

    Ω(u)={z: ||z||∞<t, ||u-z||∞≥c m(z)}

contains every E_j under consideration. Lebesgue dilation and r_j^(-n)≤m(z)^(-n) give

    ∑_{uncaptured positive j}w_j
      = α^(-n)2^(-n)∑_j∫_{E_j}r_j^(-n)dz
      ≤ α^(-n) Ψ(u),
    Ψ(u)=2^(-n)∫_{Ω(u)}m(z)^(-n)dz.                      (5.3)

The prefactor is exactly (a+1)^n/(2^n a^n): the numerator comes from the Lebesgue Jacobian and the denominator from the star normalization. This is the only replacement of each actual increment by a larger integral. It uses their pairwise disjointness and the specific constant Lebesgue background; it uses no overlap multiplicity assumption about the cubes themselves.

Combining (4.1), the captured payment, and (5.3),

    W(0)≤τ+a^(-n)+α^(-n)Ψ(u).                            (5.4)

## 6. Coordinate monotonicity puts the worst shift at a corner

We prove

    Ψ(u)≤Ψ(α,...,α) whenever ||u||∞≤α.                   (6.1)

Reflection in any coordinate leaves m, the ambient cube, and Lebesgue measure invariant, so Ψ depends only on h_i=|u_i|. Assume 0≤h_i≤α<1.

For ρ≥1 let L_i(ρ) be the length of

    [-ρ,ρ] ∩ [h_i-cρ,h_i+cρ].

Because c>1 and h_i<1≤ρ, its upper endpoint is ρ and its lower endpoint is max(-ρ,h_i-cρ). Therefore

    L_i(ρ)=min(2ρ,(1+c)ρ-h_i)
          =2ρ-[h_i-(c-1)ρ]_+.                            (6.2)

It is positive and nonincreasing as h_i increases.

Inside ||z||∞<1, m=1, and the complement of Ω is the captured intersection with the open cube u+(-c,c)^n. Its volume is ∏_iL_i(1). Thus the inner contribution to Ψ is

    1-2^(-n)∏_iL_i(1),                                  (6.3)

which is nondecreasing in every h_i.

For the outer part, parameterize the shell ||z||∞=ρ, 1<ρ<t, by the 2n coordinate faces. On the positive face in coordinate i, |h_i-ρ|<cρ always holds. The negative face is captured in that coordinate exactly when h_i+ρ<cρ, equivalently h_i<(c-1)ρ. Hence the captured face sum is

    S(ρ,h)=∑_{i=1}^n N_i(ρ,h_i)∏_{l≠i}L_l(ρ),
    N_i(ρ,h_i)=1+1_{h_i<(c-1)ρ}.                         (6.4)

For each coordinate h_j, every summand in S depends on it through either the nonincreasing nonnegative factor N_j (when i=j) or the nonincreasing positive factor L_j (when i≠j). All other factors are nonnegative. Thus S is nonincreasing in each h_j.

Here the shell decomposition has volume Jacobian 1, not Euclidean surface Jacobian. Partition points outside the unit cube according to their unique largest absolute coordinate; on one piece write z_i=±ρ and keep the remaining n−1 coordinates as independent variables in (-ρ,ρ). The change of variables is linear with absolute determinant 1. Ties between largest absolute coordinates lie in finitely many hyperplanes and have n-dimensional Lebesgue measure zero. The full shell face sum is 2n(2ρ)^(n-1). Consequently

    Ψ(u)=1-2^(-n)∏_iL_i(1)
       +2^(-n)∫₁ᵗ ρ^(-n)
          [2n(2ρ)^(n-1)-S(ρ,h)]dρ.                       (6.5)

Strict cap inequalities are used for the captured set, the complement of the ≥ comparison domain. Equalities in the remaining coordinates do not change the face lengths. The opposite-face equality h_i=(c-1)ρ affects at most one radius for each i because c>1, so it causes no ambiguity in the integral. Formula (6.5) holds also for n=1, with empty products equal to 1.

The inner term and the integrated uncovered shell term in (6.5) are both nondecreasing in every h_j. This proves (6.1). The corner is in the closed comparison box; actual centers lie in its interior, and need not attain it. No differentiation of a moving intersection volume has been used.

## 7. Evaluate the corner with the correct one-face factor

The parameter assumptions give

    c-1=(a/t-2)/(a+1)>0,
    α/(c-1)=a t/(a-2t)>t.

Therefore α>(c-1)ρ for all 1≤ρ≤t. At h_i=α each coordinate has exactly one captured face, and

    L(ρ)=(1+c)ρ-α,
    S_corner(ρ)=n[(1+c)ρ-α]^(n-1).

Also 1+c=α(2+1/t), so

    L(ρ)/(2ρ)=α g_t(ρ),
    L(1)/2=α q/2.

For 1≤ρ≤t one has 0<g_t(ρ)≤1, with equality g_t(t)=1 only at the upper endpoint. Therefore 0<α g_t(ρ)<1. Substitute the corner quantities in (6.5). The normalized full outer shell is n dρ/ρ. The captured part is

    2^(-n)ρ^(-n) n L(ρ)^(n-1)dρ
      =(n/2)[α g_t(ρ)]^(n-1)dρ/ρ.

This produces the coefficient n/2 and yields exactly

    Ψ(α,...,α)
      =1-(α q/2)^n+n log t
        -(n/2)α^(n-1)∫₁ᵗg_t(ρ)^(n-1)dρ/ρ
      =H_n(a,t).                                        (7.1)

In particular H_n(a,t)>0: its inner part 1-(α q/2)^n is positive since q/2<1, and its outer integrand is

    n[1-(α g_t(ρ))^(n-1)/2]/ρ ≥ n/(2ρ)>0.

Combining (5.4), (6.1), and (7.1), and undoing the harmless normalization,

    W(x)≤τ+a^(-n)+α^(-n)H_n(a,t)=τ+R_n(a,t)              (7.2)

for every x and every finite selected family. Both the self term and the finite-a factor α^(-n) remain present.

## 8. Recover the full operator for every L¹ input and every positive scale

Let f∈L¹ and λ>0; put E_λ={M_n f>λ}. For fixed R>0, its centered cube average is continuous in x, since its translation difference is at most

    (2R)^(-n) || |f|(·+h)-|f| ||₁ → 0.

L¹ translation continuity is valid for arbitrary L¹ functions. The supremum over all R>0 is lower semicontinuous, so E_λ is open.

Take an arbitrary compact K⊂E_λ. For each x∈K, choose an actual R_x>0 for which the average of |f| over Q(x,R_x) is >λ. Let the supplied small cube be Q(x,R_x/a); its a-star is exactly Q(x,R_x). These small open cubes cover K. Compactness provides a finite subcover, to which Sections 3–7 apply. Every selected star retains its actual bad-average property. Thus

    λ∑_j v_j
      ≤ ∑_j [v_j/|aQ_j|]∫_{aQ_j}|f(y)|dy
      = ∫ |f(y)|W(y)dy
      ≤ [τ+R_n(a,t)]||f||₁.                              (8.1)

All sums are finite and nonnegative, so there is no interchange problem. The bound (3.4), the finite subcover of K, and (8.1) give

    λ|K|≤[1+A_n(a,t)/τ][τ+R_n(a,t)]||f||₁.

Lebesgue inner regularity takes the supremum over compact K⊂E_λ and gives the same inequality for |E_λ|, even if its finiteness was not known beforehand. This proves (2.1) for arbitrary f and λ. Open or closed cubes produce identical averages because cube boundaries have Lebesgue measure zero.

This step is what promotes the finite geometric argument to the full operator. The radius R_x can be any positive number, with no upper cutoff, lower cutoff, discretization, or restriction to finitely many witness scales in the resulting operator theorem.

## 9. Take a scalar limit, then optimize the threshold

Hold n, t>1, and τ>0 fixed. The operator and every input and level set in (2.1) stay fixed while a→∞ through finite a>2t. Only the scalar right side changes. Since α→1,

    A_n(a,t)→q^n,
    R_n(a,t)=a^(-n)+α^(-n)H_n(a,t)→B_n(t).

This follows directly from the explicit formula for H_n; its integral I_n(t) is independent of a. Therefore

    C_n≤(1+q^n/τ)(τ+B_n(t)).                             (9.1)

There is no convergence requirement on subcovers, selected families, radii, increments, weights, or inputs. For every finite a the preceding theorem already applies to the same full operator.

To verify strict positivity without circularly using a square root, note q/2<1 and 0<g_t(ρ)≤1. Rewrite

    B_n(t)=1-(q/2)^n
       +n∫₁ᵗ[1-g_t(ρ)^(n-1)/2]dρ/ρ.

The inner term is positive and the integral is at least (n/2)log t>0. Hence B_n(t)>0. Put A=q^n and B=B_n(t). For every τ>0,

    (1+A/τ)(τ+B)=A+B+τ+AB/τ.

Choosing the admissible τ=sqrt(AB) yields (sqrt(A)+sqrt(B))². Substitution in (9.1) proves (1.1), including n=8.

## 10. An exact algebraic form for dimension eight

This section is optional for the proof but records a reproducible exact scalar representation. Put p=1+1/(2t), h=(1+1/t)/2. The binomial theorem and elementary integration give

    I_8(t)=∫₁ᵗ[p-1/(2ρ)]^7 dρ/ρ
      =p^7 log t+∑_{j=1}^7 binom(7,j)p^(7-j)(-1/2)^j
                                (1-t^(-j))/j,
    B_8(t)=1-h^8+8 log t-4 I_8(t).                       (10.1)

For rational t>1 all nonlogarithmic coefficients here are rational. In particular the log coefficient is 8−4p^7; no sign of that coefficient should be assumed without checking the chosen t. Formula (10.1) is an identity, not a numerical enclosure or an optimization certificate.

## 11. Audit conclusion and boundaries

The finite-a proof closes the self payment, covers rejected cubes with b+1, confines unpaid increments using b−1, handles all equal-radius contributors by the last-index tie rule, proves corner monotonicity face by face, obtains n/2 with the Lebesgue shell Jacobian, and then passes through compact level subsets to arbitrary L¹ inputs and all positive common scales. The final a→∞ operation is a limit of already proved scalar operator constants.

The hypotheses are fixed axis-parallel cubes, one common radius, and constant Lebesgue ambient density. The proof is independent of the shape or regularity of the input f beyond L¹, but its volume identities are not asserted for variable density. No non-Lebesgue theorem, lower bound, optimized numeric value, or public priority assertion is part of this result.
