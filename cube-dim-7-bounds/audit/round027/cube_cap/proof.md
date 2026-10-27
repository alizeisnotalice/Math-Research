# R027 independent cube-cap mechanism audit

**Verdict: PASS for the mechanism with the finite-a self term retained, a>2t and t>1.** I first tried to invalidate the proposed cap containment, the corner maximization, and the shell coefficient. The independent reconstruction below closes all three. No author proof or other auditor verdict was read before this proof was derived and frozen. The author's clarification confirmed only that its parameter specification is identical to the parent's. No numerical evaluation, optimization, literature search, MATLAB, commit, push, or lower-bound work was performed.

## 1. Exact statement and finite-a formula

Let n be a positive integer and let

M_n f(x)=sup_{r>0}(2r)^(-n)∫_{x+(-r,r)^n}|f(y)|dy

for arbitrary f∈L1(R^n), with Lebesgue measure. Let c_n be its least weak(1,1) constant for strict level sets. Fix t>1, a>2t and τ>0, and put

q=1+1/t,  b=a q,
α=a/(a+1),  c=(b-1)/(a+1)>1,
g_t(ρ)=1+1/(2t)-1/(2ρ),
J_n(t)=∫_1^t g_t(ρ)^(n-1) dρ/ρ.

Define

H_n(a,t)=1-α^n(q/2)^n+n log t-(n/2)α^(n-1)J_n(t),
R_n(a,t)=a^(-n)+α^(-n)H_n(a,t).

The independently proved finite bound is

c_n <= [1+(q+1/a)^n/τ][τ+R_n(a,t)].              (F)

In particular the self term a^(-n) is present before the a→∞ limit. The cap integral H_n is nonnegative because it is an actual normalized weighted uncovered volume, as derived below.

For every t>1 set

B_cap,n(t)=1-(q/2)^n+n log t-(n/2)J_n(t).

Then

c_n <= [1+q^n/τ][τ+B_cap,n(t)],                  (L)

and therefore

c_n <= [q^(n/2)+sqrt(B_cap,n(t))]^2.             (O)

Here B_cap,n(t)>0 for every t>1. No decimal estimate or optimizer is asserted. The n=7 target is exactly the specialization of these formulas, with the same common radius in all coordinates and every r>0.

## 2. Finite selection and rejected-union factor

For a finite family of open cubes Q_β=Q(y_β,r_β), choose selected cubes Q1,...,QN recursively. The increments and weights are

I_j=Q_j\(Q1∪...∪Q(j-1)),
v_j=|I_j|,
w_j=v_j/|aQ_j|=v_j/(2^n a^n r_j^n).

The I_j are pairwise disjoint measurable subsets of their small cubes Q_j. They are not increments of the star-cube union. Their measures sum to the measure of the selected small-cube union, and their nonnegative weights never change after selection.

At stage k an unselected cube centered at y is eligible if

∑_{j<k}w_j1_{bQ_j}(y)<=τ.

Choose an eligible cube of maximal radius; stop if none remains. Selected radii are nonincreasing because eligibility can only be lost when nonnegative summands are added. A finite procedure suffices, including zero increments and equal radii.

For every unselected Q_β,

∑_{j:r_j>=r_β}w_j1_{bQ_j}(y_β)>τ.               (2.1)

To prove this, consider the first selected cube whose radius is smaller than r_β, if one exists. At its selection Q_β was still available and had larger radius, so it was already ineligible, paid by previous cubes of radii>=r_β. If there is no such smaller selection, at termination all selected cubes have radii>=r_β and Q_β is ineligible. Later additions only increase the sum, proving (2.1).

For x∈Q_β and a contributing j in (2.1),

||x-y_j||_∞ <r_β+b r_j <=(b+1)r_j.

Integrating the resulting charge on the rejected union gives

|∪all Q_β| <=[1+((b+1)/a)^n/τ]∑_jv_j
=[1+(q+1/a)^n/τ]∑_jv_j.                         (2.2)

The b+1 factor is necessary and is retained exactly.

## 3. All contributors except the self term are earlier

For fixed x put W(x)=∑_jw_j1_{aQ_j}(x). If W(x)=0 there is nothing to prove. Otherwise choose k of smallest radius among its positive contributors, and among those of that smallest radius choose the **last selected index**. Translate x to0 and rescale all coordinates by r_k^(-1), so r_k=1. Weights and incidences are unchanged by this normalization.

Every positive contributor now has r_j>=1. If another positive contributor had j>k, nonincreasing radii would give r_j<=1, hence r_j=1, contradicting that k is the last positive contributor of minimum radius. Thus every other positive contributor has j<k. This tie rule is essential; choosing the first minimum-radius index would not justify the payment below.

The self term is bounded by

w_k<=|Q_k|/|aQ_k|=a^(-n).                        (3.1)

Among j≠k, call a contributor captured if y_k∈bQ_j. These captured terms are a subset of the previous terms in the eligibility test at y_k, so their total weight is <=τ. There is no self term in that eligibility sum.

For an uncaptured contributor j, y_k∉bQ_j gives

||y_k-y_j||_∞>=b r_j.                            (3.2)

Also 0∈aQ_k and 0∈aQ_j imply ||y_k||_∞<a and ||y_j||_∞<a r_j. If r_j>=t, then

||y_k-y_j||_∞<a+a r_j<=a(1+1/t)r_j=b r_j,

contradicting (3.2). Hence all uncaptured positive contributors satisfy

1<=r_j<t.                                        (3.3)

This range includes radius1 ties, but all such contributors are earlier than k by its last-index choice.

## 4. Cap containment and weighted volume

For any y∈I_j of an uncaptured contributor, I_j⊂Q_j and (3.2) imply

||y_k-y||_∞ >=||y_k-y_j||_∞-||y-y_j||_∞
>(b-1)r_j.                                       (4.1)

Set z=y/(a+1), u=y_k/(a+1). Then ||u||_∞<α. Because 0∈aQ_j and y∈Q_j,

||z||_∞<r_j<t.

With m(z)=max(1,||z||_∞), (3.3) gives r_j>=m(z); and (4.1) gives

||u-z||_∞>c r_j>=c m(z).                         (4.2)

Thus the transformed increments E_j=(a+1)^(-1)I_j lie in the common domain

Ω(u)={z:||z||_∞<t, ||u-z||_∞>=c m(z)}.

They remain pairwise disjoint. Use Lebesgue change of variables and r_j^(-n)<=m(z)^(-n) to get the pointwise-in-x uncaptured-weight estimate

∑uncaptured w_j
=((a+1)/a)^n/2^n ∑_j∫_{E_j}r_j^(-n)dz
<=α^(-n)/2^n ∫_{Ω(u)}m(z)^(-n)dz.               (4.3)

The determinant factor (a+1)^n, the star denominator a^n, and the cube-volume factor 2^n are all explicit. No density or overlap multiplicity beyond disjoint small-ball increments is assumed. Replacing the actual strict inequality in (4.2) by >= only enlarges the domain and is safe.

## 5. Corner is the worst shift: a coordinatewise proof

Define Ψ(u)=2^(-n)∫_{Ω(u)}m(z)^(-n)dz for ||u||_∞<=α. Coordinate reflections preserve cube volumes and the integrand, so the captured-volume formulas depend only on h_i=|u_i|, with0<=h_i<=α<1.

For any ρ>=1, let ℓ_i(ρ,h_i) be the length of the intersection of [-ρ,ρ] with [u_i-cρ,u_i+cρ]. Because c>1 and |u_i|<1<=ρ,

ℓ_i(ρ,h_i)=2ρ-(h_i-(c-1)ρ)_+.                   (5.1)

It is nonnegative and nonincreasing in h_i.

### Inner cube, m=1

The captured subset of (-1,1)^n is its intersection with u+(-c,c)^n and has volume

V_in(u)=∏_iℓ_i(1,h_i).

This is nonincreasing in every h_i. Therefore the uncovered inner volume is nondecreasing in every h_i and is largest at h_i=α for all i.

### Outer cube shells, m=ρ

For fixed ρ>=1, decompose the boundary of [-ρ,ρ]^n into its2n coordinate faces. On either face perpendicular to coordinate i, the remaining captured coordinates contribute area ∏_{j≠i}ℓ_j(ρ,h_j). The number of its two faces whose fixed i-coordinate satisfies the strict cap inequality |u_i-z_i|<cρ is

s_i(ρ,h_i)=1+1_{h_i<(c-1)ρ}.                     (5.2)

Indeed the face on the same-sign side as u_i is always captured; the opposite face is captured precisely when ρ+h_i<cρ. Equality excludes that opposite face for the strict captured set. Thus the captured shell-face area is exactly

S(ρ,u)=∑_{i=1}^n s_i(ρ,h_i)∏_{j≠i}ℓ_j(ρ,h_j).  (5.3)

Both s_i and ℓ_i are nonincreasing in h_i and nonnegative. For each coordinate h_l, every term of (5.3) depends on it through s_l (if i=l) or through ℓ_l (if i≠l). Hence **every term**, and consequently S, is nonincreasing in each h_l. The captured shell area is therefore smallest at a corner, and the uncovered shell area largest there. This proves the asserted corner maximization without differentiating an intersection volume as its radius changes.

The shell Jacobian is1: partition R^n outside the unit cube by its unique largest absolute coordinate, writing that coordinate as ±ρ and the remaining coordinates in (-ρ,ρ). Ties between largest coordinates form a Lebesgue-null set and do not change the volume integral. Thus the total shell area is2n(2ρ)^(n-1), and its captured part is (5.3), with no factor sqrt(n). This also proves that integrating the pointwise shell-area comparison against positive ρ^(-n)dρ gives the weighted-volume comparison. For n=1 the products over j≠i are empty products equal to1, and the same formulas hold.

We have established

Ψ(u)<=Ψ((α,...,α)).                              (5.4)

Actual normalized centers satisfy ||u||_∞<α, while the corner belongs to its closure. Extending to the closed parameter cube for an upper bound is legitimate; all formulas remain valid at h_i=α.

## 6. Exact corner volume and the one-half coefficient

The condition a>2t implies

c-1=(a/t-2)/(a+1)>0,
α/(c-1)=at/(a-2t)>t.

Consequently α>(c-1)ρ for every1<=ρ<=t. At the corner every coordinate has exactly one captured face, and

ℓ(ρ)=(c+1)ρ-α,
S_corner(ρ)=n[(c+1)ρ-α]^(n-1).                  (6.1)

The inner captured fraction is [(c+1-α)/2]^n. Normalize all volumes by2^n. The normalized full inner volume is1, the normalized full shell contribution with weightρ^(-n) is n log t, and the normalized captured shell contribution is

n/2^n ∫_1^t[(c+1)-α/ρ]^(n-1)dρ/ρ
=(n/2)∫_1^t[(c+1)/2-α/(2ρ)]^(n-1)dρ/ρ.          (6.2)

This is the source of the factor n/2: the corner captures one face per coordinate, and the remaining powers of2 are absorbed into the (n-1)st power. Counting two captured faces at the corner or using a shell Jacobian other than1 would be wrong.

Now c+1=α(2+1/t) and c+1-α=α q. Therefore (6.2) equals (n/2)α^(n-1)J_n(t), and

Ψ((α,...,α))
=1-α^n(q/2)^n+n log t-(n/2)α^(n-1)J_n(t)
=H_n(a,t).                                       (6.3)

All open/closed cap-boundary issues are accounted for by the strict face test (5.2). There is no corner-regime transition inside[1,t] for finite a>2t. At the a→∞ limit the transition reaches only the endpointρ=t, which does not change the ordinary integral.

Combining (3.1), the captured payment, (4.3), (5.4), and (6.3) gives the global weighted-overlap bound

W(x)<=τ+a^(-n)+α^(-n)H_n(a,t)=τ+R_n(a,t).        (6.4)

## 7. General L1 and all radii

For fixed r>0 the average of |f| over x+(-r,r)^n is continuous in x. Its translation increment is uniformly bounded by (2r)^(-n)|||f|(·+h)-|f|||_1, which tends to0. Hence the supremum over all r>0 is lower semicontinuous and E_λ={M_n f>λ} is open for eachλ>0.

Take any compact K⊂E_λ. For each x∈K choose an actual R_x>0 with average over Q(x,R_x) greater thanλ. The small cube Q(x,R_x/a) is open and contains x, and its star is exactly Q(x,R_x). Choose a finite subcover of K and apply the finite selection mechanism.

For the selected cubes, star averages>λ imply

λ∑_jv_j <=∫|f(y)|W(y)dy <=[τ+R_n(a,t)]||f||_1.

Using (2.2) bounds λ|K| by the right side of(F) times||f||_1. Lebesgue inner regularity of E_λ gives the same inequality for E_λ by taking the supremum over compact K. This argument also establishes finite level-set measure rather than presupposing it.

The input is arbitrary L1, with absolute value already included. The maximal operator uses allr>0; only a finite subcover for each compact level subset is selected. There is no truncation of scales or support, no atomic restriction, and no replacement of the level set by a finite witness surrogate in the final estimate. Open and closed cube averages agree because their boundaries have Lebesgue measure0.

## 8. Scalar limit and threshold optimization

Fix t>1 andτ>0 and let finite a>2t tend to∞. The operator M_n, its input, its full radius set, and its level sets stay fixed. Only the scalar bound in(F) changes, so no selector or operator limit is required.

The formula (6.3) makes the limit elementary: J_n(t) does not depend on a; α→1 and a^(-n)→0. Thus R_n(a,t)→B_cap,n(t) and(q+1/a)^n→q^n, giving(L). Moreover B_cap,n(t)>0: its defining normalized uncovered volume in the limiting corner contains the positive inner uncovered fraction1-(q/2)^n for t>1, and its shell contribution is nonnegative. Alternatively the same inner-plus-shell formula is derived directly atα=1,c=1+1/t for1<=ρ<t.

For A=q^n>0 and B=B_cap,n(t)>0,

(1+A/τ)(τ+B)=A+B+τ+AB/τ
>=(sqrt(A)+sqrt(B))^2,

with equality atτ=sqrt(AB). Since(L) holds for everyτ>0, it gives(O). This can be implemented with that fixedτ while passing a→∞; no exchange of infima is needed. Taking an infimum over t>1 is then a legitimate scalar bound. Numerical certification of any chosen t, or minimization over t, is outside this audit.

## Critical conditions and scope of this PASS

- The fixed shape is the axis-parallel cube with one common radius, and the ambient measure is Lebesgue. The face-product corner argument is cube-specific.
- Pick the last selected index among minimum-radius positive contributors; otherwise later equal-radius contributors are not paid.
- Keep the self term a^(-n) in the finite-a overlap; only its scalar limit is zero.
- Use b+1 for rejected-cube coverage and b-1 for the cap reverse-triangle containment.
- Require a>2t and t>1 before taking limits. The corner regime then follows automatically; it is not an extra assumption imposed on the inputs.
- The integral coefficient is n/2 after normalization, and shell-face ties affect no ambient volume.
- The finite-sum weights are volumes of pairwise disjoint unstarred increments. No disjointness of star cubes is used.

Under these conditions I found no break in the proposed mechanism. This is a rigorous analytic reconstruction of the mechanism, not numerical evidence and not a claim that the specialization was previously published. The author source identity, when supplied, is recorded separately after this independent audit freeze.
