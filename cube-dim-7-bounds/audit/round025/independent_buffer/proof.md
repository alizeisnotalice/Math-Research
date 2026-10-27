# R025 independent audit: free buffer and selection threshold

**Verdict: PASS, subject to the exact hypotheses and strict parameter restrictions stated below.** This is an independent reconstruction from the proposed parameters, not an inference from the author's or another auditor's acceptance. No R025 author proof was read before deriving and freezing this audit. No literature search, numerical optimization, MATLAB, or mathematical scalar computation was performed in this subtask. Numerical size and numerical minimization over t are outside this audit.

## Statement proved

Let n be a positive integer. Let D⊂R^n be an open bounded convex body symmetric about0, with nonempty interior, and let |.| denote Lebesgue measure. Its Minkowski functional is a norm and D is its open unit ball. Define

M_D f(x)=sup_{r>0} |rD|^(-1)∫_{x+rD}|f(y)|dy,

on arbitrary f∈L1(R^n). Denote by c(D) its least weak(1,1) constant, with strict level sets {M_D f>λ}. For every a>1,t>1,τ>0 set

b=a(1+1/t),
A_a(t)=((b+1)/a)^n=(1+1/t+1/a)^n,
B_a(t)=((a+1)/a)^n(1+n log t).

Then

c(D) <= [1+A_a(t)/τ][τ+B_a(t)].                     (1)

Consequently, for every t>1,τ>0,

c(D) <= [1+(1+1/t)^n/τ][τ+1+n log t],              (2)

and optimizing only the positive scalar τ gives

c(D) <= [(1+1/t)^(n/2)+sqrt(1+n log t)]^2.          (3)

One may take infimum of the right side over t>1. No value of this infimum or decimal evaluation is claimed here. In particular (3) applies to n=7 and D=(-1,1)^7, retaining one common radius and every r>0.

## 1. Finite selection, independent of f

Take a finite family of open D-balls Q_α=B_D(y_α,r_α), all radii positive. Choose balls Q1,...,QN recursively as follows. Once Q1,...,Q(k-1) have been selected, fix the **small-ball increments**

I_j=Q_j\(Q1∪...∪Q(j-1)),  v_j=|I_j|,
w_j=v_j/|aQ_j|=v_j/(|D|a^n r_j^n).

Here aQ_j is dilation about its own center y_j; the center itself is not multiplied by a. The increments are measurable and pairwise disjoint, with

∑_j v_j=|∪selected Q_j|.                            (4)

They are not increments of ∪aQ_j. Later selections never change an earlier I_j or its nonnegative weight w_j.

An unselected ball with center y is eligible when

∑_{j<k}w_j1_{bQ_j}(y) <=τ.                         (5)

Choose an eligible remaining ball of maximal radius; stop if none remains. This is a finite procedure. Eligible sets can only shrink because all added terms are nonnegative. Hence selected radii are nonincreasing. Equal radii may be ordered arbitrarily.

## 2. Rejected union charge and its exact factor

For any unselected Q_α one has

∑_{j:r_j>=r_α} w_j1_{bQ_j}(y_α)>τ.                (6)

Indeed, if a ball of radius smaller than r_α is ever selected, consider the first such selection. Q_α is still available and has larger radius, so it must be ineligible at that moment. All balls selected earlier have radius>=r_α. If no smaller ball is ever selected, all selected balls have radius>=r_α, and termination makes Q_α ineligible. This proves (6), including ties.

For any x∈Q_α and any contributing j in (6), the norm triangle inequality yields

||x-y_j||_D < r_α+b r_j <=(b+1)r_j.

Thus on the union U_rej of unselected balls,

∑_j w_j1_{(b+1)Q_j}(x)>τ.

Integration, the exact scaling |sQ_j|=s^n|Q_j|, and (4) give

|U_rej| <=(1/τ)∑_j w_j|(b+1)Q_j|
=A_a(t)/τ * ∑_jv_j.

Adding the selected union gives the finite covering estimate

|∪all Q_α| <=[1+A_a(t)/τ]∑_jv_j.                 (7)

In particular the covering dilation is b+1. Using b, or replacing b+1 by a+1, would be incorrect. No extra factor from the number of balls or from n is lost in (7).

## 3. Weighted star overlap

Define

W(x)=∑_j w_j1_{aQ_j}(x).

If W(x)=0 the desired estimate is trivial. Otherwise choose k of smallest radius among the selected balls contributing positively at x. Translate x to0 and rescale the entire finite family by r_k^(-1), so r_k=1. These operations preserve the selected order, the criterion (5), the norm-ball incidences, and the weights because numerator and denominator scale by the same power n. All positive contributions to W(0) now have r_j>=1, and

||y_k||_D<a,
||y_j||_D<a r_j for each positive contributor j.   (8)

Split W(0) into large radii r_j>=t and small radii 1<=r_j<t.

### Large range: why the self term is absent

Because t>1=r_k, every large contributor satisfies r_j>r_k. Nonincreasing radii therefore imply j<k. In particular j=k does not appear in the large range. For each such j,

||y_k-y_j||_D < a+a r_j
<=a(1+1/t)r_j=b r_j.

The eligibility condition (5) used to select Q_k thus pays every large contributor and yields

W_large(0)<=τ.                                    (9)

This is the exact place where t>1 is required. Directly setting t=1 would incorrectly place the self term in a sum controlled only before selection of k. No such endpoint substitution is made.

### Small range: cumulative disjoint increments

Let J be the finite set of positive contributors with1<=r_j<t, and define

F(u)=∑_{j∈J:r_j<=u}v_j,  1<=u<=t.

For any j counted in F(u), (8) and I_j⊂Q_j give

I_j⊂Q_j⊂B_D(0,(a+1)u).

The increments in this cumulative sum are pairwise disjoint; therefore

F(u)<=|D|(a+1)^n u^n.                             (10)

For1<=r<t the ordinary calculus identity

r^(-n)=t^(-n)+n∫_r^t u^(-n-1)du

is exact. Sum it over the finite set J with coefficients v_j. Interchanging a finite sum and an integral, without any differentiability requirement on F, gives

∑_{j∈J}v_j r_j^(-n)
=t^(-n)F(t)+n∫_1^t F(u)u^(-n-1)du
<=|D|(a+1)^n[1+n log t].                          (11)

The possible atom at radius1 is included in the t^(-n)F(t) boundary term and in the integral identity; no missing F(1) subtraction is permissible. Equality radii at interior u affect only measure-zero integration points. There is no radius-t atom in J, but including it in the cumulative upper bound would only weaken the estimate. Formula (11) is a finite sum identity, so Stieltjes regularity assumptions are unnecessary.

Divide by |D|a^n to obtain

W_small(0)<=B_a(t).

Together with (9), undoing normalization,

W(x)<=τ+B_a(t) for every x.                       (12)

## 4. General L1, all radii, actual level sets

Let f∈L1(R^n), use its nonnegative absolute value, and let λ>0. For any fixed r>0, convolution with the bounded kernel |rD|^(-1)1_{rD} is continuous: its change under an x-translation h is bounded uniformly in x by

|rD|^(-1)|||f|(·+h)-|f|||_1,

which tends to0 by L1 translation continuity. Thus M_D f, a supremum over all r>0 of continuous functions, is lower semicontinuous, and E_λ={M_D f>λ} is open.

Take any compact K⊂E_λ. At each x∈K, the strict level inequality supplies an actual R_x>0 with

∫_{x+R_xD}|f| >λ|R_xD|.

Choose the small ball Q_x=B_D(x,R_x/a). It is open and contains x; its star aQ_x is exactly x+R_xD. Compactness supplies a finite subcover of K by these small balls. Apply the finite selection argument with the fixed a,t,τ.

Every selected ball has star average>λ. Multiplication by w_j and finite summation imply

λ∑_jv_j <=∫|f(y)|W(y)dy
<=[τ+B_a(t)]||f||_1.

(The weak inequality safely includes zero increments.) Combine with (7) to get

λ|K| <=[1+A_a(t)/τ][τ+B_a(t)]||f||_1.

Lebesgue inner regularity of the open set E_λ gives the same estimate for E_λ by taking the supremum over all compact K⊂E_λ. This also establishes that its measure is finite; finiteness was not assumed in the compact approximation.

This is an inequality for every L1 input and every λ>0 for the genuine supremum over all r>0. The finite family arises only from the compact subcover. There is no input discretization, support truncation, finite scale restriction, limit of atomic measures, radial-input assumption, or level-set surrogate. Taking D=(-1,1)^7 gives the requested centered axis-parallel cube model. Closed cube boundaries have Lebesgue measure0, so closed-vs-open averaging sets produce identical averages for all L1 inputs and radii.

## 5. The a→∞ passage is scalar, not an operator limit

Fix t>1 and τ>0. Statement (1) has already been established, with the same fixed operator M_D, for every finite a>1. Its proof may choose a different finite cover and a different selector for each a; these auxiliary objects no longer occur in the resulting bound.

Let a tend to∞ in the finite scalar right side:

A_a(t)→A_∞(t)=(1+1/t)^n,
B_a(t)→B_∞(t)=1+n log t.

For any fixed f,λ the left side λ|E_λ| is independent of a. Taking this ordinary numerical limit in the universal inequality yields (2), including equality <= at the limiting constant. Equivalently c(D) is bounded by every finite-a constant and hence by their limit. No convergence or measurability of the family of selectors is required. No maximal operator, body, input, radius set, or level set is changed in this passage.

For A=A_∞(t)>0, B=B_∞(t)>0, expand

(1+A/τ)(τ+B)=A+B+τ+AB/τ.

The elementary identity

τ+AB/τ-2sqrt(AB)=(sqrt(τ)-sqrt(AB/τ))^2>=0

shows the minimum over τ>0 is attained at τ=sqrt(AB), giving (sqrt(A)+sqrt(B))^2. This is (3). One may fix precisely this τ while sending a→∞, so no exchange of two infima is needed. Finally, because the inequality holds for every t>1, taking an infimum over this scalar parameter is legitimate, whether or not the infimum is attained.

## Audit boundaries and provenance

PASS applies to the displayed analytic formulas, with n positive integer, nondegenerate symmetric convex D, finite positive a>1,t>1,τ>0 before the scalar limit. The restrictions t>1 and the exclusion of j=k in the large range are essential. The claim is an independently proved refinement of the Stein–Strömberg selection mechanism, not a claim that its numerical specialization was printed in the original1983 article. This audit establishes no decimal estimate or certified optimizer for n=7.

The original1983 article was read in the earlier R021 source task and is not re-searched here. Cached original PDF SHA-256: ecb7c3bbbf804545fd68795f5628d242ba8d03f2a3ebeab67806af9c8da6c615. This audit's argument is self-contained and does not rely on the R021/R023 constants or any other auditor verdict. Author proof identity, if received, is recorded separately for provenance and not used as a premise. No lower-bound files or previous proofs were edited.
