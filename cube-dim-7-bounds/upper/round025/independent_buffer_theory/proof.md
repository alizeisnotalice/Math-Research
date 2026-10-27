# R025: independent proof of the free-buffer selection bound

Date: 2026-10-27. Native Sol, pure theory, zero numerical calculation, MATLAB, or parameter scan. This file is authored only in the assigned new directory. Other sources, previous rounds, and the lower branch are unchanged. No commit or push. The parent retains final proof and integration responsibility; this frozen proof is supplied for another independent audit before ledger adoption.

**Conclusion.** The proposed free buffer is valid in the parameter range below. In particular the limiting constant is obtained from already proved uniform scalar inequalities, without passing a selector, cover, kernel, or operator through a limit. A simple exact corollary gives C7<2233/75<30. No approximate decimal or optimized parameter is certified here.

## 1. Precise theorem

Let n>=1 be an integer and let D be a nonempty bounded open convex origin-symmetric subset of R^n. Lebesgue measure is denoted by |.|. Put B(x,r)=x+rD and

    M_D f(x)=sup_(r>0) |rD|^(-1) integral_(x+rD) |f(y)|dy.

For any a>1, a<b<2a, and tau>0, define

    t=a/(b-a)>1,
    A_(a,b)=((b+1)/a)^n,
    B_(a,t)=(1+1/a)^n(1+n log t).

Then every f in L1(R^n) and every lambda>0 satisfy

    lambda |{M_D f>lambda}|
       <= [1+A_(a,b)/tau][tau+B_(a,t)] ||f||_1.              (1)

Equivalently, choose any a>1,t>1,tau>0 and set b=a(1+1/t). Then

    C_(n,D)<= [1+(1+1/t+1/a)^n/tau]
              [tau+(1+1/a)^n(1+n log t)].                   (2)

It follows that, for every t>1,tau>0,

    C_(n,D)<= [1+(1+1/t)^n/tau][tau+1+n log t],              (3)

and hence, for every t>1,

    C_(n,D)<= [(1+1/t)^(n/2)+sqrt(1+n log t)]^2.             (4)

The norm induced by D is used only for triangle inequalities. The operator always averages over the same centered body with one common scale; general L1 input and all positive scales are retained. In the target case D=(-1,1)^7, these are exactly centered axis-parallel cubes with one shared radius.

## 2. Finite selection with a free eligibility expansion

Start with any finite family of open D-balls. Recursively select balls B_j=B(y_j,r_j), and when selecting the jth ball define and freeze

    I_j=B_j \ union_(i<j) B_i,
    v_j=|I_j|,
    w_j=v_j/|aB_j|=v_j/(|D|a^n r_j^n).

Here aB_j means dilation about y_j. The I_j are increments of the **small** balls. They are disjoint, and sum v_j is the measure of the selected small-ball union. No star-union increments are substituted.

Among the remaining balls whose centers y satisfy

    sum_(j already selected) w_j 1_(bB_j)(y)<=tau,           (5)

select one of maximal radius. Stop when no eligible ball remains. The initial empty sum is eligible. Adding terms can only remove eligibility, so selected radii are nonincreasing. The procedure terminates on the finite supplied family. Ties and zero increments are harmless.

This constructs a family depending on a,b,tau; the arguments below use the disjointness of that actual family, without assuming its increments agree with those of any other parameter choice.

## 3. Rejection cover: the expansion is b+1

Consider any unselected ball B_alpha with radius r_alpha. Immediately before the first selected ball of radius strictly smaller than r_alpha, or at termination if there is no such selected ball, B_alpha must already be ineligible. Otherwise the maximal-radius choice could not take a smaller ball or terminate. At that stage all selected radii used are >=r_alpha. Thus

    sum_(j:r_j>=r_alpha) w_j 1_(bB_j)(y_alpha)>tau.

For x in B_alpha and a counted j,

    ||x-y_j||_D<r_alpha+b r_j<=(b+1)r_j.

Therefore sum_j w_j 1_((b+1)B_j)(x)>tau on the unselected union. Integrate and use the exact volume scaling to obtain

    |union(unselected balls)|
      <= tau^(-1) sum_j w_j |(b+1)B_j|
      = (A_(a,b)/tau) sum_j v_j.

The selected union has measure sum v_j, giving

    |union(all supplied balls)|
      <= [1+A_(a,b)/tau] sum_j v_j.                         (6)

This proves the stated cover constant, with no hidden dependence on the number of balls or on their smallest radius.

## 4. Weighted star load and the essential cutoff t>1

Fix x and put W(x)=sum_j w_j 1_(aB_j)(x). If W(x)=0 there is nothing to prove. Otherwise choose a contributing selected ball B_k with the smallest radius, translate x to 0 and normalize r_k=1. The weights are invariant under this common translation and dilation. All contributing radii are >=1, and ||y_k||_D<a.

Split W(0) into contributions r_j>=t and 1<=r_j<t. In the first group, t>1 implies r_j>r_k. The nonincreasing selected radii therefore force j<k: the self-term j=k cannot occur, and no later equal-radius term can occur. Moreover

    ||y_k-y_j||_D<a+a r_j<=b r_j,

because r_j>=t=a/(b-a). The preceding strict inequality puts y_k inside the open bB_j even when r_j=t. Hence every term in this large-radius group is contained in the selection sum (5) used at y_k, so its total is <=tau.

For the remaining contributing indices J, define

    F(s)=sum_(j in J:r_j<=s) v_j,  1<=s<=t.

If r_j<=s and 0 is in aB_j, then ||y_j||_D<a r_j, and I_j is contained in B_j. Thus all increments counted by F(s) are disjoint subsets of the ball centered at 0 of radius (a+1)s. It follows that

    F(s)<=|D|(a+1)^n s^n.                                   (7)

This is the same geometric observation for the actual new selection. The eligibility expansion b does not enter (7): the only inputs are star membership, small-ball containment, and disjoint increments.

The finite scalar identity r^(-n)=t^(-n)+n integral_r^t s^(-n-1)ds gives

    sum_(j in J) v_j r_j^(-n)
      =t^(-n)F(t)+n integral_1^t F(s)s^(-n-1)ds
      <=|D|(a+1)^n(1+n log t).

There is no continuity assumption on F and no omitted endpoint atom: r=1 is included in the finite identity, whereas r=t belongs to the large-radius group. Dividing by |D|a^n bounds the small-radius contribution by B_(a,t). We have proved pointwise

    W(x)<=tau+B_(a,t).                                      (8)

The restriction a<b<2a, equivalently t>1, ensures that the large-radius argument never charges the self-term. No assertion about t<=1 is made.

## 5. General L1 and all positive common scales

Fix f in L1, replace it by |f|, and fix lambda>0. For a fixed r the centered average is continuous in x, by L1 translation continuity and the bounded averaging kernel. Thus E={M_D f>lambda} is open.

For any compact K contained in E, choose at each x in K an actual bad radius R_x>0. Supply the small ball B_x=B(x,R_x/a), whose star aB_x is the original centered ball with average >lambda. These small open balls cover K; choose a finite subcover and apply the finite selection at the given a,b,tau. Each selected star is still bad, so

    lambda sum_j v_j
      <=sum_j w_j integral_(aB_j)|f|
      =integral |f| W
      <=[tau+B_(a,t)]||f||_1.

Together with (6), this proves (1) for K. Inner regularity gives the same inequality for E, without a prior assumption that E has finite measure. The argument uses no scale lattice, fixed radius window, finite source support, radiality, product input, or weak-type composition. In the cube application open and closed cubes have identical L1 averages.

## 6. Constant infimum as a tends to infinity

Fix n,D,t>1,tau>0. For every a>1, b=a(1+1/t) is in (a,2a), so (2) is already proved for the unchanged operator M_D and for all f,lambda. Its right-side coefficient tends to the coefficient in (3) as a tends to infinity. In fact each positive factor decreases as 1/a decreases.

For any fixed f,lambda the left side does not depend on a and is finite by (2) at one finite a. Taking the infimum of the proved scalar bounds gives (3). Equivalently one may use a sequence of finite admissible a values tending to infinity. The finite cover, number of witnesses, chosen radii, order, and increments may all change with a. No convergence or uniform cardinality for any of those objects is needed or claimed. They are eliminated before the scalar infimum is taken.

In (3), set A0=(1+1/t)^n and B0=1+n log t, both positive. Expanding its coefficient gives

    A0+B0+tau+A0 B0/tau.

It is minimized over admissible tau>0 at tau=sqrt(A0 B0). The minimum is (sqrt(A0)+sqrt(B0))^2, which proves (4). This optimization is performed after establishing (3) for each fixed tau; no parameter-dependent operator or selector limit is exchanged.

## 7. Exact rational n=7 corollary, no numerical optimization

Take t=5 and retain the fixed rational threshold tau=20/3 in (3). Then

    A0=(6/5)^7<18/5,

as the exact integer comparison 6^7=279936<281250=18*5^6 shows. The positive exponential series gives

    exp(5/3)>1+5/3+25/18+125/162+625/1944
             =10009/1944>5,

and hence log 5<5/3. Therefore B0=1+7 log 5<38/3. The coefficient in (3) is increasing in A0 and B0 for this positive fixed tau. Consequently

    C_7 <= (1+A0/(20/3))((20/3)+B0)
        < [1+(18/5)(3/20)] [20/3+38/3]
        = (77/50)(58/3)
        =2233/75<30.                                       (9)

All comparisons in this corollary are exact analytic inequalities written in the proof. No software evaluated these constants, and no numerical parameter tuning was performed. Formula (4) is the stronger accepted symbolic theorem; this branch does not certify the proposed approximate value near 29.1.

## Provenance and frozen scope

The free eligibility buffer b was proposed by the parent. The finite proof above was reconstructed independently, including all parameter restrictions, the small-radius cumulative identity, and the all-input/all-scale transfer. It uses the selection method of Stein–Strömberg already checked against the original article in R021; neither the free b theorem nor its finite constant is represented as an explicit statement printed in that article.

Only this mechanism and the exact corollary (9) were considered. The branch does not assert novelty in the public literature, optimality, or publication status. The proof is now frozen for the parent to coordinate independent review and any separately authorized scalar calculation. No automatic scan or follow-up is launched.
