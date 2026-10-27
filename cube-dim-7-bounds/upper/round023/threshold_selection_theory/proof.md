# R023: a free eligibility threshold in the finite selection mechanism

2026-10-27. Pure Sol derivation. No mathematical numerical computation, optimization, MC, MATLAB, commit or push by this theory owner. The separately authorized scalar owner may evaluate only the eight registered pairs after source locking. This proof changes the admission threshold in the Stein–Strömberg finite selection argument; the weights retain their original definition throughout.

## The parameter theorem

Let D be a bounded open convex body symmetric about0 with nonempty interior in R^n, n>=1. For all a>1, tau>0, all f inL1(R^n), and all lambda>0, its original centered homothetic maximal operator obeys

    M_D f(x)=sup_(r>0)|rD|^-1 integral_(x+rD)|f(y)|dy,

    lambda|{M_D f>lambda}|<=K_n(a,tau)||f||1,
    K_n(a,tau)=(1+A_n(a)/tau)(tau+B_n(a)),
    A_n(a)=(1+2/a)^n,
    B_n(a)=(1+1/a)^n(1+n loga).                             (1)

The supremum defining M_D ranges over every common radius r>0. For fixed a, minimizing the proved expression over tau gives

    min_(tau>0)K_n(a,tau)=(sqrt(A_n(a))+sqrt(B_n(a)))²,
    tau_opt=sqrt(A_n(a)B_n(a)).                              (2)

For the7Dcentered axis-parallel cube, a64 and the rational tau6 already give the completely analytic bound

    C7<1040858/21375<49.                                    (3)

This is a direct weak(1,1)upper bound, not a conditional Green certificate. No assertion of optimality or literature novelty is made. Attribution of the underlying selection mechanism to Stein–Strömberg1983and its primary section2is retained in sources.json. Their fixed threshold1does not constrain the independently proved positive threshold here.

## 1. Finite construction: weights are not rescaled

Use the norm induced by D. Given a finite family of balls B(y,r)=y+rD, set B*=B(y,ar), B**=B(y,(a+1)r), B***=B(y,(a+2)r), all with the same center. Select recursively a remaining ball of largest radius among those whose center y satisfies

    sum_(j<k) w_j1_(B_j**)(y)<=tau.                          (4)

After selecting B_k put

    I_k=B_k\union_(j<k)B_j,
    w_k=|I_k|/|B_k*|.                                       (5)

These are the increments of the SMALL-ball union. They are disjoint, and sum|I_j|=|union_selected B_j|. Their volumes and weights stay fixed after definition. Adding later nonnegative weights only removes eligibility, so the selected radii are nonincreasing: a remaining ball eligible later was eligible earlier, and the earlier choice maximized radius. The first selection is possible for every nonempty family because its empty sum0is <=tau. Termination after finitely many choices is automatic; zero-volume increments cause no problem.

No step requires tau<=1, an integer tau, or tau bounded below by an individual weight. The estimates below are uniform over all finite families even though changing tau changes the selected subfamily and its increments.

## 2. Rejection charge and the first factor

For an omitted ball B_alpha with radius rho and center y_alpha,

    sum_(j:r_j>=rho)w_j1_(B_j**)(y_alpha)>tau.                (6)

To prove this, if a smaller radius is ever selected, look just before the first such choice. All earlier selected radii are >=rho. Since B_alpha was not chosen, it was ineligible; otherwise it would have a strictly larger eligible radius than the chosen ball. If there is no smaller selected radius, termination says B_alpha is ineligible with all selected radii >=rho. Equal-radius ties and later nonnegative contributions do not change(6).

For x in B_alpha, each contributing j in(6) has

    ||x-y_j||_D <rho+(a+1)r_j <=(a+2)r_j.

Thus V(x)=sum_j w_j1_(B_j***)(x)>tau on the union of omitted balls. Integrating this pointwise inequality yields

    |union_omitted B_alpha|<=tau^-1 integral V
      =A_n(a)/tau *sum_j|I_j|.

Adding the selected union gives

    |union_all B_alpha|<=(1+A_n(a)/tau)sum_j|I_j|.           (7)

The denominator tau comes precisely from rejection. Replacing the weight by an unmentioned scaled quantity in this step would alter its meaning.

## 3. Weighted star overlap: large part tau, small part independent of tau

Let W(x)=sum_jw_j1_(B_j*)(x). If it is nonzero, choose a contributing selected ball B_k with minimal contributing radius. Translate x to0 and dilate by r_k, so r_k=1, ||y_k||_D<a, and every contributing radius is >=1. The weights are unchanged by this common dilation.

For a contributing j with r_j>=a, strict a>1 implies j<k. Moreover

    ||y_k-y_j||_D <a+a r_j<=(a+1)r_j,

where r_j>=a is essential. Those terms occur in(4)when k was selected, so their total is <=tau.

The remaining contributing set J has1<=r_j<a. For1<=t<=a define F(t)=sum_(j in J,r_j<=t)|I_j|. Every counted increment lies in B(0,(a+1)t), because its parent ball contains0after a-expansion. Disjointness alone gives

    F(t)<=|D|(a+1)^n t^n.                                  (8)

Importantly, this estimate uses only containment and disjointness. It is true for whichever subfamily the tau-dependent algorithm selected and contains no tau. No selection-density assumption has been hidden in its right side.

The finite identity r^-n=a^-n+n integral_r^a t^(-n-1)dt gives, including any atom at radius1,

    sum_(j in J)|I_j|r_j^-n
      =a^-nF(a)+n integral_1^aF(t)t^(-n-1)dt
      <=|D|(a+1)^n(1+n loga).

Dividing by |D|a^n proves the small contribution is <=B_n(a). Therefore

    W(x)<=tau+B_n(a) for every x.                           (9)

This argument does not need disjoint selected balls, disjoint star balls, differentiability of F, or a continuous radius distribution. The exact finite identity fixes the boundary sign and retains the radius1atom.

## 4. Original all-radius centered maximal operator on generalL1

Take h=|f|. For each fixed r>0 its centered average is continuous in x, since its translated difference is bounded uniformly by |rD|^-1||tau_z h-h||1. Thus E={M_D f>lambda} is open. If K subset E is compact, for each x in K choose an actual witness radius R_x>0 with integral_(x+R_xD)h>lambda|R_xD|. The smaller open balls x+(R_x/a)D contain their centers and cover K, so a finite subcover exists. No measurable selection of radii is required.

Apply the preceding construction to that finite cover. For every selected ball its star is a witness with average>lambda. Multiplying by its actual small increment volume and summing gives

    lambda sum_j|I_j|<=integral_Rn h(y)W(y)dy
                      <=[tau+B_n(a)]||f||1.

The integral over all space is finite by the pointwise bound(9)and h inL1. This handles unbounded support and integrable tails without any truncation. Combining with(7)gives(1)for K. Supremizing over compact subsets by Lebesgue inner regularity gives(1)for E, even if its measure was not known finite beforehand. The constant is independent of the selected finite family, its largest or smallest radius, K, f and lambda. Hence the argument includes every common positive radius and every generalL1input.

Closed cube averages coincide with open cube averages because their boundaries are Lebesgue null. A weak-type convention with{M>=lambda}follows by thresholds increasing tolambda and continuity of measure after obtaining finite upper bounds. No comparison with an auxiliary maximal operator incurs an extra multiplicative factor; the construction directly uses centered witnesses for the original operator.

## 5. Why threshold normalization does not remove the gain or degenerate it

The exact product is

    K_n(a,tau)=A+B+tau+AB/tau.                              (10)

Scaling weights to w_tilde=w/tau would make the eligibility threshold1, but the actual witness charge would then read lambda sum|I_j|<=tau integral h sum w_tilde1_(B_j*). Its small-overlap bound would be B/tau, not B. Multiplying all factors gives exactly(10)again. The genuine freedom is the choice of selected subfamily admitted by(4), not a false algebraic gain from changing units.

Since A,B>0, the terms tau and AB/tau diverge at the respective ends tau infinity and tau0. AM-GM gives tau+AB/tau>=2sqrt(AB), with equality precisely at the positive tau in(2). The theorem is uniform in f,lambda,so choosing this tau as a function only of n,a is legitimate. There is no dependence on unknown increments, function data, time, spatial point or radii in this optimization.

## 6. Independent purely rational bound below49

For n7,a64 put A=(33/32)^7 and B=(65/64)^7(1+42log2). The elementary binomial/geometric estimate(1+x)^7<=(1-7x)^-1 for0<=x<1/7 gives

    A<=32/25, (65/64)^7<=64/57.

The positive exponential series gives exp(7/10)>12013/6000>2, so log2<7/10 and

    B<(64/57)(152/5)=9728/285.

At the rational threshold tau6, positivity permits

    K_7(64,6)<(1+(32/25)/6)(6+9728/285)
                =(91/75)(11438/285)
                =1040858/21375<49,

where49*21375=1047375>1040858. This proves(3)without evaluating a logarithm or running a scalar program. The upcoming finite numerical batch can sharpen the displayed rational upper bound, but is not needed to establish the improvement from the free threshold.

## 7. Finite next calculation and the limit of the claim

The registered eight pairs are a in{32,36,40,64}, tau in{6,7}, n7. Each is a valid bound independently; a rigorous minimum of their outward upper endpoints remains a valid bound. These are a fixed finite sample of the proved parameter family, not a claim of global optimization in a or tau. No repeats of the oldtau1baseline or new large Green work are needed. Exact square-root optimization in tau is proved analytically, but evaluating additional optimized rows is not part of the eight-pair batch.

This is a direct upper estimate for the same centered cube weak-type constant. It does not determine that best constant, establish a lower bound, or prove novelty against all literature. The primary source's fixed-threshold mechanism and the new free-threshold derivation are distinguished. Parent and independent audit retain responsibility for adoption and numerical reporting.
