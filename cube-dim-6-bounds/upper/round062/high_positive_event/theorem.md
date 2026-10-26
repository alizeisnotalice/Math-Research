# High positive event: a full-time two-event budget

This is a proof prepared before production. No new scalar has been run.
Keep the original body32 equal-time Cauchy model, q=6/5 and rho=1/64.
On F={(1,y):y in[0,1]^5}, let H be any subset of the union of the31
unsigned tubes T_A={dist(x,span(1_A,1_B))<rho}, with A containing1,
1<=|A|<=5 and mean_B(x)>3/4. Cap, low and middle exclusions may be
retained in H or safely relaxed. Fix the new parameter-event threshold
b=1/4; it is not the old cap threshold1/2.

## 1. High-coordinate and all remaining-subspace distances

In a tube write x_A=alpha1_A+e and x_B=u1_B+z with zero-sum residuals.
Then |e|^2+|z|^2<rho^2, x1=1 gives alpha=1-e1, and xi<=1 gives
e1>=ei. B coordinates are >3/4-rho=47/64. A coordinates are
xi=1+ei-e1>=1-sqrt(2)rho>47/64. These statements include singleton
blocks and coordinates on the face boundary. Therefore every x in H
has every coordinate >=47/64.

For any |a|<b, let y=x-a ones. Every coordinate of y is then

    yi>=47/64-b=d=31/64>0.                              (1)

The following uniform inequalities hold for all such y:

    |y|^2>=6d^2,
    dist(y,span v)^2>=(10/3)d^2 for every body v!=ones,
    dist(y,span(v,w))^2>=2d^2
                   for distinct body v,w neither ones. (2)

For the line bound, write y=d ones+delta with delta>=0. If v has h
minus signs, 1<=h<=5, its squared distance at d ones is
(2/3)h(6-h)d^2. The linear increment has nonnegative coefficients
2d[1-v_i(6-2h)/6], and the remaining quadratic term is a squared
projection norm. The minimum integer value h(6-h) is5.

For the pair bound, split coordinates into the agreement and
disagreement supports of v,w. Both are nonempty; the plane has one
signed support direction on each. If both sign patterns were constant,
that plane would contain ones. A plane spanned by two distinct body
lines contains exactly those two body lines: a candidate body vector
av+bw must satisfy |a+b|=|a-b|=1, hence ab=0. Thus a plane containing
ones would have ones as one of its two original body lines. That is
impossible after omitting ones. At least one support block consequently
has opposite signs. For two such coordinates the least squared
residual is (yi+yj)^2/2>=2d^2, proving the last bound in(2).

This is an analytic all-line/all-pair argument, not a new enumeration.
It does not claim that two remaining signed lines are geometrically
disjoint or that the original high tubes are disjoint.

## 2. Exact positive event split and large-parameter branch

Let r_t be the density of the sum of the31 parameters other than ones,
and f_t(a)=t/[pi(a^2+t^2)]. Split the actual density as

    p_t=P_t+B_t,
    P_t(x)=integral_(|a|>=b) f_t(a)r_t(x-a ones) da,
    B_t(x)=integral_(|a|<b) f_t(a)r_t(x-a ones) da.        (3)

Both are nonnegative, the events are disjoint and exhaustive up to a
null threshold boundary, and every other parameter is unrestricted.
There is no assumption that all other jumps are small. The all-t
pointwise coarea proof in accepted A041/F041 is unchanged by choosing
b=1/4: the restored quotient along ones has the same exact Q5 peak,
and f_t(a)<=t/[pi(b^2+t^2)] on the first event. Thus

    P_t<=D t^(-4)/(b^2+t^2),
    D=R1/(32pi^6), R1=67320683/3251298631680.            (4)

There is one selected parameter, hence no32-fold deletion factor.
The pointwise coarea equality, including its uniform integrable line
tail argument, is read from the accepted pointwise_coarea_note; it is
not inferred merely from continuity of the density.

## 3. Full complementary event, with the correct extension order

Inside the31-parameter density, apply the accepted A041 positive
deletion inequality three times. At quotient dimension k=0,1,2, the
remaining frame operator is 32I-(P ones)(P ones)^T<=32I and there are
31-k columns. This produces the three time coefficients992,960,928.
On the restricted event |a|<b, use(2) to pull the denominators outside
the a-integral. Only after this replacement extend the still-positive
integral, retaining its factor f_t(a), to all real a. It then restores
the original033 three-dimensional quotient density and its accepted
Q3 peak. The extension is valid even if ones lies in the deleted triple
span; then its projection is zero and integral f_t=1.

The already accepted complete sum of26970 ordered triples excluding
ones gives Crem=44743079/(2970pi^6). No new quotient peak or old triple
enumeration is used. The complete bound for every x in H and t>0 is

    B_t(x)<=Crem/[(6d^2+992t^2)
                   ((10/3)d^2+960t^2)(2d^2+928t^2)].   (5)

The positive-coordinate gap is required only before denominators are
replaced. It is not incorrectly asserted for the added |a|>=b points.
All actual parameter events, including multiple large parameters and
cancellations, are already included in(3)-(5).

## 4. True high-union volume and the P short-time estimate

For k=|A|, m=6-k, use the genuine face coordinates of034/054:
dy=sqrt(km) du de dz, normal dimension4, and maximal-A-coordinate cone
fraction1/k for the radial volume integrand. Retain 3/4<u<=1, a length
1/4, and remove other box clipping only after taking nonnegative
integrals. This gives

    |H_A|<=(pi^2/8)rho^4 sqrt(m/k).

The all31 positive volume union, allowing any overlaps, therefore gives

    |H|<=V:=(pi^2/8)rho^4
                    sum_(k=1)^5 binom(5,k-1)sqrt((6-k)/k). (6)

This is an upper volume formula, not the exact volume of clipped H.
Singleton A or B blocks are included. Later norm estimates use this
one whole-union V, not31 copies of the same large-parameter budget.

The first coordinate of the full body32 law is Cauchy(32t), and the full
density is even in each of the other five coordinates. Its positive
face slice is t/[pi(1+1024t^2)]. Since P_t<=p_t and H is a subset of
the positive face,

    integral_H P_t<=t/[pi(1+1024t^2)]<=t/pi.             (7)

No such reflection symmetry is assumed for P_t itself. Put C=D/b^2.
Interpolation between(7) and P_t<=C t^(-4), followed by Minkowski, gives

    ||integral_0^epsilon P_t dt||_q
        <=S(epsilon):=(6/7)C^(1/6)pi^(-5/6)epsilon^(7/6). (8)

## 5. Full P tail, one analytic epsilon, and full B time integral

For general b>0, a direct partial fraction identity gives

    t^(-4)/(b^2+t^2)
      =b^(-2)t^(-4)-b^(-4)t^(-2)+b^(-4)/(b^2+t^2).

Therefore its full tail is

    T_b(epsilon)=1/(3b^2 epsilon^3)-1/(b^4 epsilon)
                          +atan(b/epsilon)/b^5.        (9)

The derivative is the negative integrand and the limit at infinity is
zero. The positive arctangent branch is used. Consequently

    ||integral_epsilon^infinity P_t dt||_q
       <=V^(5/6)D T_b(epsilon).                        (10)

Use only the fixed analytic choice

    epsilon=((D/b^2)pi V)^(1/5).                       (11)

It uniquely minimizes the loose expression
S(epsilon)+V^(5/6)C/(3epsilon^3): its derivative changes sign once at
epsilon^(25/6)=(C pi V)^(5/6). This does not claim to optimize the exact
tail(9), scan parameters or change the underlying equal-time model.

For the complement put c=(sqrt6 d,sqrt(10/3)d,sqrt2 d),
A=(992,960,928), and alpha_i=c_i/sqrt(A_i)>0. The elementary identity

    integral_0^infinity product_i(t^2+alpha_i^2)^(-1) dt
      =pi(alpha0+alpha1+alpha2)/
        [2alpha0 alpha1 alpha2(alpha0+alpha1)
                    (alpha0+alpha2)(alpha1+alpha2)]

follows from the two-factor integral pi/[2ac(a+c)] and partial
fractions, and extends to coincident positive roots by dominated
convergence. Pulling out product A_i gives J_B, the full integral of
the denominator reciprocal in(5). Thus

    ||integral_0^infinity B_t dt||_q<=V^(5/6)Crem J_B.   (12)

## 6. Complete new high budget and proof status

The two parameter events are disjoint and exhaustive; the two P time
pieces are disjoint and exhaustive, and B is integrated once on the
entire positive axis. Their spatial densities need not be disjoint.
Minkowski for finite positive time truncations and monotone convergence
therefore prove the prepared formula

    ||G||_Lq(H)<=S(epsilon)
       +V^(5/6)D T_b(epsilon)+V^(5/6)Crem J_B.          (13)

This new event bound does not use the large ones-first Phi coefficient
from060 and does not rerun the042 cap scalar. It applies even if the
whole high union still contains its cap intersections, because(1)
holds there too; subtracting prior spatial cells only improves it.

All constants and physical choices are fixed in registration.json.
Actual production requires a separate explicit parent release and
hash binding. Each producer/checker must be killed and reaped at a
180-second timeout, with no automatic rerun. No scalar result or final
PASS is claimed at this preparation stage. The parent retains the
mathematical proof and integration responsibility. Middle long times,
far, global g, whole-face I and C6 remain outside this result.
