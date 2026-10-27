# R023 independent audit: free selection threshold

Date: 2026-10-27. Native Sol; pure mathematical audit, zero numerical or MATLAB runs, no commit or push. This file belongs only to `audit/round023/threshold_selection`. Earlier sources and the lower branch are untouched. The parent retains final proof responsibility and ledger integration.

**Verdict: accept the free-threshold theorem below.** This audit re-proves the finite selection, its two bounds, and the operator transfer. It does not infer correctness merely by replacing a symbol in R021 or by citing a companion auditor. The dependence of the selected family and increments on tau creates no missing factor.

## Full theorem

Let n>=1 be an integer and D a nonempty bounded open convex origin-symmetric subset of R^n. For f in L1(R^n) define the centered maximal operator

    M_D f(x)=sup_(r>0) |rD|^(-1) integral_(x+rD) |f(y)|dy.

For every a>1 and tau>0 put

    A=(1+2/a)^n,
    B=(1+1/a)^n(1+n log a).

Then, for every lambda>0,

    lambda |{M_D f>lambda}|
       <= (1+A/tau)(tau+B) ||f||_1.                         (T1)

Consequently the optimal weak-type constant satisfies

    C_(n,D) <= (sqrt(A)+sqrt(B))^2                            (T2)

for each a>1, by the legitimate choice tau=sqrt(AB). This statement includes general signed or complex L1 inputs, all positive common dilation scales, and arbitrary D of the stated class. For D=(-1,1)^7 it applies to the exact target operator. No rectangular maximal operator or independent coordinate scales are introduced.

## 1. Actual selection for arbitrary positive tau

Take an arbitrary finite supplied family of open D-balls. If B_1,...,B_(k-1) have been selected, write

    I_j=B_j \ union_(i<j) B_i,
    v_j=|I_j|,
    w_j=v_j/|aB_j|,

where aB_j is the dilation about its own center. Among remaining balls whose centers y satisfy

    sum_(j<k) w_j 1_((a+1)B_j)(y)<=tau,                     (S)

choose one of maximal radius. Stop when none is eligible. The empty initial sum is eligible because tau>0. The weights and I_j are fixed at their selection times. They are not recomputed against a later union, and they are not divided by tau.

Adding a selected term only increases every remaining eligibility sum. Therefore eligibility can only be lost, and the selected radii are nonincreasing. This holds for every positive tau, not only tau=1. For every realized selection the I_j are pairwise disjoint, with sum v_j equal to the measure of the selected small-ball union. Although the family and increments depend on tau, these facts are unconditional properties of the constructed order.

## 2. Covering estimate: the factor is exactly 1+A/tau

Let B_alpha be unselected with radius r_alpha. Just before the first selected ball with smaller radius, B_alpha cannot be eligible: otherwise choosing that smaller ball violates the maximal-radius rule. If there is no such smaller selected ball, the same conclusion holds at termination. Thus, using only selected radii at least r_alpha,

    sum_(j:r_j>=r_alpha) w_j 1_((a+1)B_j)(y_alpha)>tau.

Equal radii and zero increments cause no exception. For x in B_alpha, the triangle inequality implies that every contributing term also has x in (a+2)B_j. Consequently

    sum_j w_j 1_((a+2)B_j)(x)>tau

on the unselected union. Integration gives

    |union(unselected balls)|
       <= tau^(-1) sum_j w_j |(a+2)B_j|
       = (A/tau) sum_j v_j.

Adding the selected union gives

    |union(all supplied balls)| <= (1+A/tau) sum_j v_j.      (C)

The only inverse tau enters here, from integrating the rejection threshold. The selected union contributes the separate 1.

## 3. Star load: the factor is exactly tau+B

Fix x and define W(x)=sum_j w_j 1_(aB_j)(x). If W(x)=0 the desired estimate is immediate. Otherwise choose a contributing selected ball B_k of smallest radius, translate x to 0, and scale r_k to 1. The weights are unchanged by this normalization. All contributing radii are >=1 and ||y_k||_D<a.

For contributing radii r_j>=a, one has j<k and

    ||y_k-y_j||_D<a+a r_j<=(a+1)r_j.

Thus this part of W(0) is a sub-sum of the left side of (S) at y_k and is at most tau. Positivity justifies discarding the other eligibility terms. The cutoff r_j>=a is essential; using r_j>=1 instead would be invalid.

For the remaining contributing indices J, with 1<=r_j<a, define

    F(t)=sum_(j in J:r_j<=t) v_j,  1<=t<=a.

Every such I_j lies in B_j and its center satisfies ||y_j||_D<a r_j. Hence the increments counted by F(t) are disjoint subsets of the D-ball centered at 0 of radius (a+1)t. Therefore

    F(t)<=|D|(a+1)^n t^n.                                   (V)

Neither (V) nor disjointness uses the selection threshold. The exact finite-sum identity

    sum_(j in J) v_j r_j^(-n)
      =a^(-n)F(a)+n integral_1^a F(t)t^(-n-1)dt

follows by summing r_j^(-n)=a^(-n)+n integral_(r_j)^a t^(-n-1)dt. It includes any atom at r=1. The endpoint r=a is assigned to the large-radius part. Applying (V) and dividing by |D|a^n bounds the small-radius part by B. No tau appears, and no regularity assumption on the step function F is needed.

We have proved, for the family selected using that same tau,

    W(x)<=tau+B  for every x.                               (W)

In particular the proof does not insert a factor tau into the small-radius estimate, and does not require the increments to agree with those produced at tau=1.

## 4. General L1, all scales, and constant optimization

Fix f,lambda and put E={M_D f>lambda}. For each fixed radius the average is a continuous function of the center by L1 translation continuity, so E is open. For any compact K contained in E and every x in K, select an actual bad radius R_x>0 and supply the small ball B_x=B(x,R_x/a). Its star aB_x is the original centered common-radius bad ball. The small balls form an open cover of K; extract a finite subcover and use the tau-selection above.

Every selected star still has average greater than lambda. Thus

    lambda sum_j v_j
      <= sum_j w_j integral_(aB_j)|f|
      = integral |f| W
      <= (tau+B)||f||_1.

Combining this with (C) proves (T1) for K. Inner regularity gives (T1) for E. This does not assume E has finite measure beforehand. No radius window, support cutoff, discrete input, or product assumption has been used. Open and closed cubes have equal L1 averages because their boundaries have zero Lebesgue measure. The bad-average threshold remains lambda throughout; tau changes only the finite geometric selection criterion.

The exact product is

    (1+A/tau)(tau+B)=A+B+tau+AB/tau.

Since A,B>0, the last two terms are minimized at tau=sqrt(AB), with minimum 2sqrt(AB). This is an attained admissible parameter choice and proves (T2), not an exchange of an unjustified limit or two unrelated selected families.

## 5. A rational seven-dimensional corollary without numerical optimization

For n=7,a=64 the quantities are

    A=(33/32)^7,
    B=(65/64)^7(1+42 log 2).

The binomial/geometric comparison (1+x)^7<=1/(1-7x) yields

    A<=32/25,
    (65/64)^7<=64/57.

The strict bound log 2<7/10 follows from

    exp(7/10)>1+7/10+49/200+343/6000=12013/6000>2.

Therefore

    B<(64/57)(152/5)=512/15.

Take the simple admissible rational threshold tau=20/3. The product in (T1) is increasing in each positive A,B with this tau fixed, so

    C_7 <= (1+A/(20/3))((20/3)+B)
        < [1+(32/25)(3/20)] [20/3+512/15]
        = (149/125)(204/5)
        =30396/625<49<128.                                  (R)

This optional corollary was obtained entirely by exact symbolic inequalities and rational arithmetic in the proof. No numerical run or search was made. The optimized expression in (T2) is the accepted sharper formula; this audit does not certify an approximate decimal or an optimized choice of a.

## Provenance and boundary of the verdict

The baseline geometric selection is the independently checked R021 adaptation of Stein–Strömberg §2. R023's arbitrary tau and optimized formula are new explicit deductions, not statements quoted as printed in the 1983 article. The parent supplied the proposed tau formula; the proof above was reconstructed independently before consulting a companion's verdict. Source integrity receipts for the original article and R021 baseline remain in this clone's `audit/round021/explicit_public_upper/`.

The verdict accepts the theorem and the rational corollary (R), without claiming optimality, novelty in the literature, or publication status. No ledger was edited. The finite invariants, source freezing, and full operator transfer are submitted to the parent for final integration.
