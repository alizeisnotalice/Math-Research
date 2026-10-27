# R021 independent audit: explicit parameter in the Stein–Strömberg selection proof

2026-10-27. Native Sol independent proof audit. No mathematical numerical run, optimization, MC, MATLAB, commit or push. I read the original Stein–Strömberg paper, section2 printed pages260–263, and the literature agent's parameter draft, but did not adopt another geometric auditor's reasoning. The proof below re-establishes each finite estimate and the generalL1/all-radius passage. Parent retains final proof responsibility.

## Result and source status

For a bounded open convex body D symmetric about0 with nonempty interior in R^n, let

    M_D f(x)=sup_(r>0) |rD|^-1 integral_(x+rD)|f(y)|dy.

For every integer n>=1 and a>1 the following bound holds:

    |{M_D f>lambda}| <= C_n(a)||f||1/lambda,
    C_n(a)=[1+(1+2/a)^n]
              *[1+(1+1/a)^n(1+n log a)].                    (1)

In particular the independent purely analytic specialization below gives

    C_7 <= C_7(64) <10013/125<128.                           (2)

Here C_7 on the left means the weak(1,1)constant for the original Lebesgue centered axis-parallel cube maximal operator on generalL1 and all positive common radii. No Green potential, source floor, moment budget, empirical point or changed operator enters this proof.

The primary paper is E.M.Stein and J.-O.Strömberg, Behavior of maximal functions in R^n for large n, Arkiv för Matematik21(1983),259–269, DOI10.1007/BF02384314. Its Theorem1 states an unspecified absolute multiple of n log n. It does not state(1) or the numeral in(2). The selection method is theirs; the free a and exact finite radius-sum estimate here give an explicit quantitative consequence. The original reproduction and local PDF/text hashes are in sources.json. This is attribution of the mechanism, not a claim that the published article printed this finite-dimensional numeral.

## 1. A finite greedy construction with immutable weights

Use the norm induced by D. For a ball B=B(y,r)=y+rD, let B*=B(y,ar), B**=B(y,(a+1)r), B***=B(y,(a+2)r); these dilations keep the center fixed. Write v_D=|D|>0.

Start with any nonempty finite family of such balls of positive radii. At stage k, with earlier selected balls B1,...,B_(k-1), choose a remaining eligible ball of largest radius. Eligibility for a ball of center y means

    sum_(j<k) w_j 1_(B_j**)(y)<=1.

After selecting B_k define its small-ball increment

    I_k=B_k \ union_(j<k)B_j, w_k=|I_k|/|B_k*|.

If there is no eligible remaining ball, terminate. All subsequent weights are nonnegative, and previously defined increments and weights never change. Thus eligibility can only be lost. In particular selected radii r1,...,rN are nonincreasing: a ball eligible later was eligible earlier, and the earlier rule chose the largest eligible radius. The finite procedure ends. Zero increments are allowed and have weight0.

The I_j are pairwise disjoint Borel sets, and

    sum_j|I_j|=|union_j B_j|.                               (3)

Only these small-ball increments are used. Replacing them by increments of B* would not be this proof.

## 2. Every omitted ball is charged to sufficiently large selected balls

Fix an unselected ball B_alpha of radius rho and center y_alpha. Then

    sum_(j:r_j>=rho) w_j 1_(B_j**)(y_alpha)>1.               (4)

Indeed if the algorithm ever selects a radius strictly smaller than rho, consider its first such stage. All earlier selected radii are >=rho. The omitted ball was still present; if its current eligibility sum had been <=1, it would outrank the chosen smaller radius. Hence the sum over those earlier large balls was >1. Future contributions cannot decrease it. If no selected radius is smaller than rho, termination makes every omitted ball ineligible with all selected radii >=rho. Ties at rho cause no problem. This proves(4) without assuming a particular tie-breaking order.

For x in B_alpha, every contribution in(4) has

    ||x-y_j||_D <rho+(a+1)r_j <=(a+2)r_j.

Hence the finite function V(x)=sum_j w_j1_(B_j***)(x) exceeds1 on the union of all omitted balls. Integrating and using volumes of dilates gives

    |union_omitted B_alpha|
      <=integral_Rn V
       =((a+2)/a)^n sum_j|I_j|.

Together with(3), this proves the first constant

    |union_all B_alpha|<=[1+((a+2)/a)^n]sum_j|I_j|.          (5)

No overlap multiplicity bound for the enlarged balls has been used here; the weights and the integral of their indicators do the charging.

## 3. Uniform control of the weighted star overlap

Let W(x)=sum_j w_j1_(B_j*)(x). If it is zero the desired estimate is immediate. Otherwise choose a contributing selected ball B_k with smallest radius among the contributing balls. Translate x to0 and rescale by r_k, so r_k=1, ||y_k||_D<a and all contributing r_j>=1. The weights are invariant under this common change of coordinates.

Split W(0)=W_large+W_small, where the first sum has r_j>=a and the second1<=r_j<a. A contributing j in the first sum occurs earlier than k, since a>1 and selected radii do not increase. Also

    ||y_k-y_j||_D <a+a r_j <=(a+1)r_j,

where the last inequality uses r_j>=a, not merely r_j>=1. The eligibility condition used when B_k was selected therefore gives

    W_large<=1.                                             (6)

This radius threshold is essential. Using the entire range r_j>=1 here would leave an invalid triangle inequality and can lose the argument's main constant.

For the small range let J be its finite contributing index set, write v_j=|I_j|, and define

    F(t)=sum_(j in J, r_j<=t)v_j, 1<=t<=a.

For a term counted in F(t), its center obeys ||y_j||_D<a r_j and I_j is contained in B_j. Thus every point of I_j has norm<(a+1)r_j<=(a+1)t. Disjointness of the I_j gives

    F(t)<=v_D(a+1)^n t^n.                                  (7)

The integral estimate can be audited without any Stieltjes convention. For every individual r in[1,a),

    r^-n=a^-n+n integral_r^a t^(-n-1)dt.

Multiply by v_j and add over the finite J. Interchanging a finite sum with the nonnegative integral gives the exact identity

    sum_(j in J)v_j r_j^-n
       =a^-n F(a)+n integral_1^a F(t)t^(-n-1)dt.

The atom at r=1 is included in both its appropriate boundary contribution and the full interval integral by this identity; no F(1)term has been discarded or counted an extra time. Using(7) now gives

    sum_(j in J)v_j r_j^-n
       <=v_D(a+1)^n[1+n log a].

Division by v_D a^n is exactly the definition of w_j on the normalized balls. Therefore

    W_small<=((a+1)/a)^n[1+n log a].

Combining with(6), for every original spatial point x,

    W(x)<=A_n(a):=1+((a+1)/a)^n[1+n log a].                 (8)

The estimate is pointwise for a finite family. It assumes neither disjoint B_j nor disjoint B_j*, only disjoint increments I_j.

## 4. All radii and arbitrary L1 input

Let f be any real or complex L1 function, use h=|f|, and fix lambda>0. For fixed r the centered average is continuous in x. For example its change under translation by z is uniformly bounded by |rD|^-1 ||tau_z h-h||1, which tends to0 by L1 translation continuity. Therefore E={M_D f>lambda}, a supremum-of-continuous-functions strict level set, is open and measurable.

Fix any compact K subset E. For each x in K there is an actual positive radius R_x with

    integral_(B(x,R_x)) h >lambda |R_xD|.

Use the smaller open ball B(x,R_x/a). These balls contain their centers and cover K, so compactness provides a finite subcover. No measurable choice of R_x is required. For every ball in this finite family, its star is precisely the original witness B(x,R_x), whose average exceedslambda.

Apply the finite lemma and retain its selected balls/increments/weights. For each selected positive increment,

    lambda|I_j| <w_j integral_(B_j*)h;

for a zero increment the nonstrict inequality holds. Summing, using the finite overlap function and then(8),

    lambda sum_j|I_j|
       <=integral_Rn h(y)W(y)dy
       <=A_n(a)||f||1.

This global integral is finite by(8) and h inL1. There is no hidden compact-support assumption or omitted tail of f. Combining with(5) and K inside the original finite cover gives

    lambda|K| <=[1+((a+2)/a)^n]A_n(a)||f||1.

Take the supremum over all compact K subset E. Lebesgue inner regularity yields(1), even if E was initially allowed infinite measure; the uniform finite right side proves it is finite. No truncation of the range of radii, smallest admissible radius, upper radius bound or support cutoff has been imposed. Each compact argument simply uses its finite family of arbitrary positive witness radii, and the constants do not depend on that family.

This direct witness construction works for the original centered maximal operator. The original article introduces an auxiliary enlarged-ball operator; its optional comparison back to M does not require an extra multiplicative constant in our argument because we never use that comparison. In particular no factor e or2is missing from(1).

If the weak-type convention uses{M_D f>=lambda}, apply(1)at thresholds increasing tolambda. The relevant sets eventually have finite measure, and continuity of measure from above gives the same bound. For D=(-1,1)^7, replacing open cubes by their closed versions leaves every L1 average unchanged, since cube boundaries have Lebesgue measure0. Thus this is the requested centered-cube operator on allL1, all common positive radii.

## 5. A completely analytic a64 specialization

For n7,a64,

    C_7(64)=[1+(33/32)^7]
                  [1+(65/64)^7(1+42 log2)].                 (9)

For0<x<1/7, binomial coefficients bounded by7^k give(1+x)^7<=sum_(k>=0)(7x)^k=(1-7x)^-1. Consequently

    1+(33/32)^7<=57/25, (65/64)^7<=64/57.

Moreover exp(7/10)>1+7/10+(7/10)²/2+(7/10)³/6
    =12013/6000>2,

so log2<7/10. Applying these upper estimates to(9), preserving positivity,

    C_7(64)<(57/25)[1+(64/57)(152/5)]
             =10013/125<128.                               (10)

All comparisons here are exact elementary inequalities. No numerical logarithm, interval arithmetic, optimization or scientific scalar batch was run by this audit owner. An even simpler log2<1 estimate would already give2809/25<128, but(10)is stronger and uses only a three-term exponential lower polynomial plus1. The result remains an explicit consequence of an established selection mechanism, not a claim that a new numerical experiment found it.

## 6. Audit conclusion and limitations

The candidate formula(1)passes this independent proof audit. Its two factors come from different finite statements: omitted-ball measure charging and weighted star overlap. The logarithm has a positive sign after exact finite integration, and the denominator of every weight is the star volume a^n|B_j|. Ties, zero increments, omitted balls, strict center eligibility and the smallest contributing radius were all treated above.

If adopted after parent/other independent review, (10)is a direct weak(1,1)upper bound for the original centered-cube problem. It requires none of the unresolved globalg assumptions in the Green route. This audit does not assert optimality, novelty relative to all literature, a lower bound, or any value of the exact best constant. Source attribution and direct proof are essential: quoting only an unspecified O(n logn)theorem would not by itself justify the finite numeral.
