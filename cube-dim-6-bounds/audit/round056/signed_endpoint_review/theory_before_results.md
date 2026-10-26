# F056 independent theory, written before opening author arithmetic results

Scope: the original body32 frame, the same 31 E_A, rho=1/64,
L=1/4, q=6/5. The F054 geometry and full-time C037 envelope are
frozen inputs. No defining M integral is to be rerun.

For b=x_B>=0 and a mixed sign vector with h<=m-h negative entries,
let P and N be the nonnegative sums on the two sign groups. Then
|P-N|<=max(P,N), and both P^2 and N^2 are at most (m-h)||b||^2
by Cauchy-Schwarz. Hence

    dist(b,span sigma)^2 >= (h/m)||b||^2 >= h u^2.

The zero-sum A residual only adds a nonnegative term to the full plane
distance. Unsigned projection dominance on b>=0 independently gives
dist(x,P_sigma)>=r. Thus dist>=max(r,sqrt(h)u). Zeros and ties cause
no exception. This estimate is applied on the physical positive domain
before enlarging the integration domain. The replacing integrand is
nonnegative and radial in the four normal coordinates, so the F054
Jacobian sqrt(km) and maximal-coordinate cone fraction 1/k give an
upper integral sqrt(m/k) 2pi^2 times the r,u integral. Positivity of
coordinates is not asserted on the enlarged ball.

For c=sqrt(h), rho<cL. Split the u integral at r/c:

    integral_0^L max(r,cu)^(-q) du
      = [q r^(1-q)/c - L^(1-q)/c^q]/(q-1).

After multiplying by r^3 and integrating r from 0 to rho, the two
terms are q rho^(5-q)/[(q-1)(5-q)c] and
-L^(1-q)rho^4/[4(q-1)c^q]. Dividing by the old unsigned radial-mean
upper formula L rho^(4-q)/(4-q) gives

    R(alpha)=q(4-q)/[(q-1)(5-q)] alpha
             -(4-q)/[4(q-1)] alpha^q,
    alpha=rho/(cL).

Thus at q=6/5, R=(84/19)alpha-(7/2)alpha^(6/5).
For 0<alpha<=1, R>=(35/38)alpha>0;
R'>=84/19-21/5=21/95>0, and R(1)=35/38<1.
These are ratios of enlarged upper formulas, not ratios of physical
clipped-domain integrals or lower estimates of G.

Modulo global negation, the mixed-plane multiplicity is C(m,h) for
h<m/2 and C(m,m/2)/2 at the even tie. The one unsigned plane plus
these classes totals 2^(m-1). For m=5,4,3,2,1 the lists are
(5,10),(4,3),(3),(1),(). Each incident plane contains exactly two
body representatives, hence two ordered pairs of class min(k,m),
by the already audited F054 analytic classification.

Let S_k and B_k be the full F054 scalar intervals for the singular
and background norm formulas. Linearity before Minkowski and
equal per-plane coefficients yield

    F_m=(1+sum_h N_mh R_h^(1/q))/2^(m-1),
    U_new,k=F_m S_k+B_k.

For m>=2, 0<F_m<1. All background intervals are copied exactly.
For m=1 the complete k5 case is copied, with no reevaluation.
The positive union inequality uses the same G and weights
(1,5,10,10,5); overlaps are allowed. A strict comparison of new
upper and the saved old054 formula lower endpoint certifies budget
improvement only. High endpoints, middle long-time tails, outside
tubes, global g/I and any C6 conclusion remain outside scope.
