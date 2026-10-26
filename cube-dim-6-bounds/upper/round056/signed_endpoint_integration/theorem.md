# Positive-B signed endpoint integration

Keep the original body32 frame, positive face F={(1,y):y in[0,1]^5},
rho=1/64, q=6/5 and L=1/4. For A containing1, k=|A|, m=6-k,
use the same E_A as054, with u=mean_B(x) in[0,L] and
r=dist(x,U_A)<rho. All time t>0 is already integrated in the frozen
C037 envelope. This task changes only the spatial estimate of its
incident signed-plane terms. The054 background is kept unchanged.

Use054 only after its final independent F054 PASS is bound to this
package. Its complete incident and background scalar intervals are
inputs; no defining M integral or complete054 budget is recomputed.

## 1. A pointwise distance bound using positivity

For an incident plane P_sigma=span(1_A,sigma_B), choose the global sign
of sigma so that it has h minus signs, 1<=h<=m/2. Let b=x_B>=0 and
write S_plus=sum_(sigma_i=+1)b_i, S_minus=sum_(sigma_i=-1)b_i.
Since these sums are nonnegative,

    |sigma.b|=|S_plus-S_minus|<=max(S_plus,S_minus)
                 <=sqrt(m-h)||b||.

The last inequality is the ordinary Cauchy-Schwarz inequality on each
of the two sign groups, whose sizes are m-h and h<=m-h. Therefore

    dist(b,span sigma)^2
       =||b||^2-(sigma.b)^2/m >=(h/m)||b||^2 >=h u^2,     (1)

where ||b||^2>=m u^2. This proof includes zero coordinates and ties.
On the full six-dimensional plane, its A residual contributes a
nonnegative additional squared distance, so dist(x,P_sigma)>=sqrt(h)u.
The unsigned projection dominance from054 also gives

    dist(x,P_sigma)>=dist(x,U_A)=r.

Consequently, on the actual nonnegative domain only,

    dist(x,P_sigma)^(-q)<=max(r,sqrt(h)u)^(-q).            (2)

The value at r=u=0 can be taken as+infinity; this null set does not
affect the finite integral. There is no claim that different mixed
planes or their neighborhoods are disjoint.

## 2. Genuine normal measure and legitimate enlargement

The054/F034 normal coordinates are x_A=alpha1_A+e,
x_B=u1_B+z, with zero-sum e,z, e1>=e_i and x1=1. Their Jacobian is
sqrt(km), and the four-dimensional normal radius is r=|(e,z)|.
For radial integrands, the A maximal-coordinate cone has fraction1/k,
including k1 with zero-dimensional A residual. At m1 no mixed plane
exists, so that type will be copied from054 without reevaluation.

Inequality(2) is first applied on the true E_A, where b>=0. Only after
replacing the original distance integrand by the nonnegative function
max(r,sqrt(h)u)^(-q) do we remove the remaining box/positivity clipping.
The replacement function is defined on the enlarged ball even though
(1) need not hold outside the original domain. Set inclusion and
nonnegativity justify the enlarged integral. Thus we do not assume
clipping is symmetric, nor apply the orthant estimate to negative b.

For c=sqrt(h), an upper bound for the true mixed-plane singular integral
is the radial enlarged-domain quantity

    J_h,k = sqrt(m/k)*2pi^2
            integral_0^rho r^3 integral_0^L max(r,cu)^(-q) du dr.
                                                                    (3)

The old enlarged unsigned quantity is

    J_0,k = L sqrt(m/k)*[2pi^2/(4-q)]rho^(4-q).            (4)

It is also the054 J upper formula. Neither expression is asserted to
be the exact integral over the clipped physical domain.

## 3. Exact one-dimensional and radial integration

For 0<r<=rho<cL and 1<q<4, split at u=r/c. Direct integration gives

    integral_0^L max(r,cu)^(-q)du
      =q r^(1-q)/[(q-1)c]-L^(1-q)/[(q-1)c^q].            (5)

Here rho=1/64<L=1/4 and h>=1, so the split always lies inside[0,L].
Although (5) is written as a difference, it is a positive integral.
Inserting into(3) gives

    J_h,k=sqrt(m/k)*2pi^2 *
        {q rho^(5-q)/[(q-1)c(5-q)]
            -L^(1-q)rho^4/[4(q-1)c^q]}.                 (6)

For q=6/5, define alpha_h=rho/(sqrt(h)L). The exact ratio of the two
enlarged upper formulas is independent of k:

    R_h=J_h,k/J_0,k
       =(84/19)alpha_h-(7/2)alpha_h^(6/5), h=1,2.        (7)

Only h1 and h2 occur because m<=5. This ratio is not claimed to equal
the ratio of actual clipped-domain integrals.

For 0<alpha<=1, positivity follows from alpha^(6/5)<=alpha and
84/19-7/2=35/38>0. Its derivative satisfies

    R'(alpha)=84/19-(21/5)alpha^(1/5)>=21/95>0.

Hence 0<R(alpha)<=R(1)=35/38<1. This proves in advance that the
proposed replacement improves each mixed-plane enlarged integral,
without needing a decimal comparison or optimizing rho.

## 4. Complete mixed-sign multiplicities

Signs sigma and -sigma define the same incident plane. After choosing
the smaller minus group, for h<m/2 there are exactly binom(m,h) classes:
each class has a unique such representative. If m is even and h=m/2,
global negation pairs distinct sign vectors, so there are binom(m,h)/2
classes. The two uniform sign choices represent one unsigned plane.
Thus

    N_mh=binom(m,h) for h<m/2,
    N_m,m/2=binom(m,m/2)/2,
    1+sum_h N_mh=2^(m-1).                               (8)

The lists (N_m1,N_m2) for m5,4,3,2,1 are
(5,10),(4,3),(3,0),(1,0),(0,0). Each plane has exactly two ordered
pairs of the same class j=min(k,m), as proved in054. Thus this count
retains all 2^m ordered incident terms, without direction enumeration.

## 5. Norm estimate using unchanged054 coefficients

Let S_k and B_k denote the exact054 upper formulas for the singular
and background Lq norm contributions. Every incident ordered term has
coefficient M_j/(r0,k^2 d0,k^2), independent of sigma. Minkowski, the
one unsigned plane, all mixed multiplicities and(3)-(8) give

    new incident norm <= F_m S_k,
    F_m=[1+sum_h N_mh R_h^(5/6)]/2^(m-1),
    U_k,new=F_m S_k+B_k.                                (9)

The use of Minkowski does not require mutual disjointness. It is an
upper bound on the sum of all positive terms in the genuine C037
Green envelope. The background retains every old nonincident term,
with exactly the complete054 saved interval and no recomputation.
For m1, F_1=1 and the entire k5 result is copied from054 unchanged.
For m>=2, (7)-(8) imply 0<F_m<1 and strict improvement of the upper
formula, without asserting a lower bound or an asymptotic value of G.

Finally use the common function G and the same31-domain positive union:

    ||G||_Lq(union_A E_A)
       <=[sum_(k=1)^5 binom(5,k-1) U_k,new^(6/5)]^(5/6). (10)

All overlaps are allowed; the exact weights remain(1,5,10,10,5).
This sharpens the same low-endpoint full-time upper budget only.
High endpoints, middle long times, outside-tube regions, the global
Green floor and a full I/g or C6 conclusion remain outside this result.
