# A selected large Cauchy parameter with a fully paid complement

Work with the original body32 frame V={(1,s2,...,s6)}, its independent
Cauchy(t) parameters a_v, and v0=(1,...,1). All parameters below have the
same time t. Define

    F={(1,y): y in[0,1]^5},
    C_rho={x in F: dist(x,span(v0))<rho},
    rho=1/64, epsilon=1/1000, a*=1/2, q=6/5.

This note proves a finite-time upper bound for the actual short potential
on C_rho. It does not estimate the rest of the face or the time tail
epsilon to infinity, and does not certify I, g or C6.

## 1. Exact positive partition by one real parameter

Let r_t be the six-dimensional density of the sum over V minus{v0}.
Since its frame Gram is32I-v0v0^T>=26I, its characteristic function is
integrable and r_t is bounded and continuous. Put f_t(a)=t/[pi(a^2+t^2)].
The actual full density splits exactly as p_t=P_t+B_t, where

    P_t(x)=integral_{|a|>=1/2} f_t(a) r_t(x-a v0) da,
    B_t(x)=integral_{|a|< 1/2} f_t(a) r_t(x-a v0) da.       (1)

Both are nonnegative continuous densities of subprobability measures.
No assertion is made that the other31 parameters are small on the first
event. Multiple large parameters, large parameters along other lines,
and cancellation events remain in (1). The two events are disjoint and
exhaust the parameter space up to a null threshold boundary.

On |a|>=1/2, f_t(a)<=t/[pi(1/4+t^2)]. Orthogonal coarea along v0 gives

    integral_R r_t(x-a v0) da
       =|v0|^-1 q_(v0perp),t(P_(v0perp)x).

The quotient law is the original body32 quotient: the omitted v0 has
zero projection. Its exact033 peak is Q5=sqrt(6)R1/(2pi)^5, with
R1=67320683/3251298631680. Thus, for all x and t>0,

    P_t(x)<= [R1/(32pi^6)] t^-4/(1/4+t^2).                 (2)

Only one parameter was removed. There is no32-term deletion sum in(2).
In particular the coefficient after dropping 1/4+t^2 to1/4 is
C_big=R1/(8pi^6). Its ratio to the old whole-face coefficient32R1/pi^6
is exactly1/256. Formula(2) bounds only the first branch, not all p_t.

## 2. Geometry for the complementary event

Write x=v0+h, h1=0, and e=P_(v0perp)h. Then h=e-e1 v0 and
|h|^2=|e|^2+6e1^2<=6|e|^2, because e1^2<=5|e|^2/6.
Consequently |h|<=eta:=sqrt(6)rho on the closed cap as well.

For |a|<=1/2 and any subspace U,

    dist(x-a v0,U)
       >=|1-a|dist(v0,U)-|h|
       >=(1/2)dist(v0,U)-eta.                             (3)

The exact geometry of the original frame gives

    |v0|=sqrt(6),
    min_(v!=v0) dist(v0,span(v))=sqrt(10/3),
    min_(v,w distinct, neither v0) dist(v0,span(v,w))=sqrt(2).

The last two claims are checked with exact Gram ratios for all31 lines
and465 unordered pairs, saved in results.json. They can also be seen
from033: the triple Gram determinants are64,96,128,160,192; the64 class
has pair determinants20 or32, and the other classes have determinant
at least96 while every pair determinant is at most36. Therefore every
relevant ratio det Gram(v0,v,w)/det Gram(v,w) is at least2. Two distinct
single-coordinate flips of v0 attain2. This uses genuine distances to
each individual plane, not the distance to the union S2.

Set the strictly positive constants

    c0=sqrt(6)/2-eta,
    c1=sqrt(10/3)/2-eta,
    c2=1/sqrt(2)-eta.                                     (4)

They bound the three distances in(3) for U={0}, one other body line,
and a pair of other body lines, respectively. In particular the cap
contains S1 and intersects its incident S2 planes, but x-a v0 stays
uniformly separated from every plane used in the next paragraph.

## 3. Three positive deletions inside the31-parameter remainder

Fix U spanned by k=0,1,2 already deleted directions, none equal to v0,
and W=Uperp. The remaining N=31-k projected columns w_j satisfy

    sum_j w_j w_j^T =32I_W-(P_Wv0)(P_Wv0)^T <=32I_W.

All are nonzero, since any three distinct body directions are independent.
For z=sum_j w_j b_j, positivity yields the actual product inequality

    product_j f_t(b_j)
      <= [32t/(pi(|z|^2+32Nt^2))]
          sum_i product_(j!=i) f_t(b_j).                  (5)

Indeed division by the positive product reduces(5) to
|z|^2+32Nt^2<=32sum_j(b_j^2+t^2). Integrating against nonnegative test
functions and changing variables along the deleted direction gives the
coarea factor1/|P_Wv_i|. Repeating three times multiplies these factors
to1/J3, J3=sqrt(det Gram(v,w,u)); there is no extra coordinate Jacobian.

The statement is pointwise, not conditioning on a zero-probability fiber.
The test-function proof first gives an a.e. density inequality. All its
quotient densities are continuous: the remaining projected frame after
omitting v0 has Gram32I_W-(P_Wv0)(P_Wv0)^T>=26I_W. Thus Fourier inversion
and dominated convergence extend the inequality to every z. This applies
at each stage before the next deletion.

Apply(5) three times to r_t(x-a v0). On |a|<1/2 and x in C_rho, (3)-(4)
bound its denominators below by

    (c0^2+992t^2)(c1^2+960t^2)(c2^2+928t^2).               (6)

KEEP f_t(a) and integrate over |a|<1/2. Once(6) has been pulled outside,
positivity permits extending that integral to all real a. For each
deleted triple U this integral is exactly

    integral_R f_t(a) q^(omit v0)_(Uperp),t(P_Uperp x-a P_Uperp v0) da
       =q_(Uperp),t(P_Uperp x).

This restores the original033 quotient law. If v0 happens to belong to
the triple span, its projection is zero and the same equality holds
because integral f_t=1. No new quotient peak, independence of the large
events, or unequal time parameter is introduced.

After applying the exact Q3 peaks, all t^3 factors cancel. We obtain

    B_t(x)<= C_rem /
      [(c0^2+992t^2)(c1^2+960t^2)(c2^2+928t^2)],          (7)

    C_rem=(32/pi)^3 sum_(ordered v,w,u in V\{v0}, distinct)
                    Q3(span(v,w,u))/J3(v,w,u).

Every non-v0 ordered triple is present. The five unordered counts are
435,1740,435,580,1305 in033 class order4110,3210,2220,3111,2211; they sum
to4495. One exact finite-frame check also verifies these counts directly.
Transitivity provides a second combinatorial derivation: each direction
belongs to3n_j/32 triples of class j, so omitting v0 leaves29n_j/32.
Since Q3/J3=R3/(2pi)^3,

    C_rem=[24576/pi^6] sum_j n'_j R3,j = (29/32)B3*.        (8)

This use of B3* is an exactly restricted sum with uniformly positive
complementary-event distances. It is not the old global common-distance
bound on the cap, whose d1 and d2 vanish there. In particular

    B_t(x)<=M_rem:=C_rem/(c0^2 c1^2 c2^2).                 (9)

Equations(2) and(7), added together, are a true pointwise upper bound for
p_t on this cap, for every t>0, including the corner and all S2 crossings.
The new scalar evaluation is restricted to0<t<=epsilon.

## 4. A finite short-potential cap budget

Let G_epsilon(x)=integral_0^epsilon p_t(x)dt. On the five-dimensional face
we use the actual coordinate measure dy. The exact034 cap volume is

    Vcap=(4sqrt(6)pi^2/45)rho^5.                           (10)

For the first branch, P_t<=p_t and the genuine positive-face slice bound
from B036 gives

    integral_Crho P_t dy <=t/[pi(1+1024t^2)]<=t/pi.

No1/32 symmetry claim is made for P_t itself; the already established
whole-density slice bound is used by domination. Interpolation with(2)
and 1/4+t^2>=1/4 gives

    ||P_t||_Lq(Crho) <= C_big^(1/6) pi^(-5/6) t^(1/6).

For the other branch use(9) and(10). Minkowski for nonnegative finite
truncations, followed by monotone convergence, proves

    ||G_epsilon||_Lq(Crho)
       <= (6/7) C_big^(1/6)pi^(-5/6)epsilon^(7/6)
          +M_rem Vcap^(5/6)epsilon.                       (11)

Thus possible infinite Green values on singular zero-measure sets cause
no gap in the Lq assertion. Both events, all other large parameters, and
the entire small-time interval are paid by(11). Its numerical evaluation
uses just rho=1/64, epsilon=1/1000, a*=1/2, with no parameter scan.

The unresolved quantities are the Green contribution at times t>epsilon
on the same cap, the S2 middle portions outside this cap, and the rest of
the face, as well as the separate global g floor and Hnorm goal. Formula
(11) alone is not a bound for the complete face integral I. No old failed
two-envelope cap budget is rerun here.
