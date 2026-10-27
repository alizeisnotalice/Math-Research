# R021 independent audit: explicit Stein–Strömberg upper bound

Date: 2026-10-27. Native Sol, pure mathematical audit; zero numerical calculations, MATLAB runs, commits, or pushes. This audit is restricted to the single frozen proof and original article §2. It does not modify their sources. The parent retains final proof responsibility and ledger integration.

**Verdict: accept the mathematical derivation.** No missing constant, invalid selection step, endpoint loss, or general-L1/all-scale gap was found. The displayed finite-dimensional formula is a quantitative derivation from the published proof, not an explicitly stated numerical conclusion of the original article. No optimality or publication-novelty claim is made.

## Accepted theorem, with its full range

Let n>=1 be an integer and let D be a nonempty bounded open convex set in R^n, symmetric about 0. Write |.| for n-dimensional Lebesgue measure, B(x,r)=x+rD for r>0, and

    M_D f(x)=sup_(r>0) |rD|^(-1) integral_(x+rD) |f(y)|dy.

For **every** a>1, every f in L1(R^n), and every lambda>0,

    lambda |{M_D f>lambda}| <= A_n(a) B_n(a) ||f||_1,

where

    A_n(a)=1+(1+2/a)^n,
    B_n(a)=1+(1+1/a)^n(1+n log a).                            (T)

The body D and input f need not be radial or products. The supremum includes every positive common dilation scale of D. Taking D=(-1,1)^7 gives the requested centered axis-parallel cubes with the same radius in all seven coordinates. Closed cubes give the same integrals for L1 input.

## Independent verification of the finite geometric lemma

Start with any finite family of open D-balls. For selected balls B_j=B(y_j,r_j), define

    I_j=B_j \ union_(i<j) B_i,
    v_j=|I_j|,
    w_j=v_j/|a B_j|=v_j/(|D| a^n r_j^n).

The dilation aB_j is about y_j. The I_j are increments of the **small** balls, and are pairwise disjoint. The weights are fixed when their balls are selected. Select a maximal-radius remaining ball whose center y obeys

    sum_(j already selected) w_j 1_((a+1)B_j)(y)<=1.

Stop when none is eligible. All added terms are nonnegative, so eligibility can only be lost. Consequently selected radii are nonincreasing. Zero increments and equal radii do not affect this argument.

For an unselected ball B_alpha with radius r_alpha, consider the stage just before the first selected radius strictly smaller than r_alpha, or termination if no such stage exists. Had its eligibility sum still been <=1, the maximal-radius rule could not select a smaller ball or terminate. Hence

    sum_(j:r_j>=r_alpha) w_j 1_((a+1)B_j)(y_alpha)>1.

If x is in B_alpha, each term counted here also has x in (a+2)B_j, by the norm triangle inequality. Integrating over the unselected union gives

    |union(all balls)| <= A_n(a) sum_j v_j.                  (G1)

The factor A_n contains the selected union itself, whose measure is exactly sum v_j. There is no further union-volume or overlap constant.

Next fix a point x with positive star load W(x)=sum_j w_j 1_(aB_j)(x). Among balls that contribute choose B_k of smallest radius. Translate x to 0 and normalize r_k to 1. Every contributing radius is now >=1 and ||y_k||_D<a. Split W(0) into radii r_j>=a and 1<=r_j<a.

For the first group, r_j>1 guarantees j<k. Also

    ||y_k-y_j||_D<a+a r_j<=(a+1)r_j,

because r_j>=a. Thus every such term is counted by the eligibility sum at y_k, so this part of W is <=1. The condition r_j>=a, rather than merely >=1, is indispensable and is satisfied by the frozen draft.

For the remaining contributing indices J define

    F(t)=sum_(j in J:r_j<=t) v_j,  1<=t<=a.

Their increments lie in the ball centered at 0 of radius (a+1)t, since I_j is contained in B_j and ||y_j||_D<a r_j. Disjointness therefore proves

    F(t)<=|D|(a+1)^n t^n.                                   (G2)

There is no need to assume the balls themselves disjoint. This cumulative bound concerns the disjoint small-ball increments, not star increments.

For each r_j in [1,a), the exact scalar identity

    r_j^(-n)=a^(-n)+n integral_(r_j)^a t^(-n-1)dt

gives, on summing the finitely many nonnegative terms,

    sum_(j in J) v_j r_j^(-n)
      =a^(-n)F(a)+n integral_1^a F(t)t^(-n-1)dt
      <=|D|(a+1)^n(1+n log a).

The atom at the smallest possible radius 1 is included in both the finite identity and F(t); the endpoint a is already assigned to the large-radius group. There is no missing Stieltjes boundary term. Dividing by |D|a^n and adding the large group proves

    W(x)<=B_n(a)  for every x.                              (G3)

The scaling normalization does not change the weights. The proof uses only finite sums and ordinary integration; no limit over radii or integration-by-parts regularity assumption is hidden in (G3).

## Transfer from the finite lemma to the actual operator

For fixed r, x -> average_(x+rD)|f| is continuous. Indeed its difference under a translation h is bounded uniformly in x by |rD|^(-1) || |f|(.+h)-|f| ||_1, which tends to zero. Hence E={M_D f>lambda} is open.

For any compact K contained in E, at each x in K choose an actual bad radius R_x>0. Supply the small ball B_x=B(x,R_x/a). Its star aB_x is the original genuine bad ball, and the B_x form an open cover of K. Extract a finite subcover and apply the finite selection above. Every selected star has average greater than lambda, so

    lambda sum v_j
      <=sum_j w_j integral_(aB_j)|f|
      =integral |f| W
      <=B_n(a)||f||_1.

Combining this with K contained in the supplied union and (G1) gives lambda |K|<=A_n(a)B_n(a)||f||_1. Inner regularity gives (T) for E, even if its measure was initially unknown to be finite.

The compact extraction restricts only the number of witnesses used to prove one compact-set inequality. It does not impose a scale lattice, a fixed radius window, finite source support, smoothness, or any extra property on f. The reduction to |f| handles signed or complex input. No composition of weak-(1,1) operators occurs.

The original article uses an auxiliary uncentered-in-the-small-ball operator. That extra operator is unnecessary here: the preceding construction supplies centered bad stars directly. Consequently no comparison factor such as e is missing from (T).

## Exact seven-dimensional specialization, checked analytically

For n=7,a=64, (T) gives

    C_7 <= [1+(33/32)^7] [1+(65/64)^7(1+42 log 2)].          (S)

For x>=0 with 7x<1, the binomial coefficients satisfy binom(7,k)<=7^k, so (1+x)^7<=1/(1-7x). Thus

    1+(33/32)^7<=57/25,
    (65/64)^7<=64/57.

The positive exponential series proves

    exp(7/10)>1+7/10+49/200+343/6000
              =12013/6000>2,

so log 2<7/10. Substitution in (S), with no decimal evaluation, gives

    C_7 <= C_7(64)
        < (57/25)[1+(64/57)(152/5)]
        =10013/125
        <81<128.                                           (S1)

All comparisons here are direct exact algebraic inequalities, not numerical evidence. The weaker draft certificate using only log 2<1 is also correct: C_7(64)<2809/25<128. Formula (S), rather than the coarse rational ceiling in (S1), is the sharper accepted analytic bound. Parameter optimization was not run and is outside this audit.

## Source correspondence and integrity

Original article: E. M. Stein and J.-O. Strömberg, *Behavior of maximal functions in R^n for large n*, Ark. Mat.21 (1983), 259-269, DOI 10.1007/BF02384314. Original reproduction: https://archive.ymsc.tsinghua.edu.cn/pacm_download/116/7398-11512_2006_Article_BF02384314.pdf . Independently read the full §2 geometric argument, printed pp.260-263, through the original PDF web extraction and the local section extraction. The original takes star dilation n, proves the selection and covering argument, and estimates the small-radius part by geometric bins. The free a and cumulative finite-integral estimate above are the audited explicit derivation.

Frozen source PDF:

    /Users/zhengzhihao/Documents/Codex/2026-10-27/n-n-5-1-1-verify/work/cube7-literature-sources/stein-stromberg1983.pdf
    SHA256 ecb7c3bbbf804545fd68795f5628d242ba8d03f2a3ebeab67806af9c8da6c615

Frozen author draft:

    /Users/zhengzhihao/Desktop/T/Math-Research-cube7-literature/cube-dim-7-bounds/upper/round021/explicit_public_upper/PROOF-DRAFT.md
    SHA256 0ab07aafedc89191fb28f632601b0be72f355cf3745ee22c4820bd865c96f042

Both hashes were checked against the actual files. No source was edited. Initial unavailable local PDF-extractor attempts produced no mathematical result; the original-source web extraction and existing local extraction were then used. This is an independent mathematical acceptance of the frozen formula, submitted to the parent for final integration, not an assertion that the original article explicitly states the finite bound 10013/125.
