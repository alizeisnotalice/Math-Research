# All31 unsigned middle domains: one common short-time interface

The model remains the original body32 equal-weight Cauchy frame. Fix
rho=1/64, epsilon=1/1000. For every nonempty proper A containing coordinate1,
let k=|A|, B=A complement, m=6-k, v=ones, w=1_A-1_B and U_A=span(v,w).
Define

    D_A={x=(1,y), y in[0,1]^5: dist(x,U_A)<rho,
                                      u=mean_B(x) in[1/4,3/4]}.

Only the short potential S_epsilon=integral_0^epsilon p_s ds is estimated.
All time parameters in the probability model are the same s.

## 1. True block geometry, including one-point blocks

Write x_A=alpha1_A+e, x_B=u1_B+z, with e,z zero-sum on their blocks.
Because x1=1, e1=1-alpha. The center is x_A(u)=1_A+u1_B and

    h=x-x_A(u)=(e-e1 1_A,z),
    |h|^2=dist(x,U_A)^2+k e1^2<=k dist(x,U_A)^2,

since e1^2<=(k-1)|e|^2/k. For k=1, e=0 and the same bound becomes
|h|<=rho. Thus eta_k=sqrt(k)rho is valid in every case.

For a vector with value p>0 on A and value q on B, the largest absolute
inner product with a body direction is kp+m|q|. It is attained by ones
if q>=0 and by w if q<0 (up to the same unoriented line). The exact
nearest-line squared distance is therefore

    kp^2+mq^2-(kp+m|q|)^2/6=(km/6)(p-|q|)^2.              (1)

This accounts for every signed body line, not only the two named ones.

On |a_v|<1/2, the center residual has p=1-a_v>=1/2 and q=u-a_v.
If q>=0, p-|q|=1-u>=1/4. If q<0, p-|q|=1+u-2a_v>=u>=1/4.
On |a_w|<1/16 it has p=1-a_w>=15/16, q=u+a_w>=3/16 and
p-q=1-u-2a_w>=1/8. Subtract eta_k using the1-Lipschitz property of
distance. Uniform first-two-distance bounds are

    h_v0=sqrt(k)(1/2-rho),
    h_v1=(1/4)sqrt(km/6)-sqrt(k)rho,
    h_w0=sqrt(k)(15/16-rho),
    h_w1=(1/8)sqrt(km/6)-sqrt(k)rho.                      (2)

All are strictly positive. In particular h_w1/sqrt(k)>0 follows from
m/384>rho^2=1/4096 for every integer m>=1; the h_v1 gate is weaker.
The other two gates are immediate rational inequalities. No case is
rescued by changing rho or the parameter thresholds.

## 2. Every other signed pair plane stays separated at original x

At x_A(u), a different unsigned partition has some block containing both
a coordinate1 and a coordinateu. Its distance to its constant-on-block
line is at least(1-u)/sqrt(2)>=1/(4sqrt(2)). For any signed pair plane
with that partition, the unsigned projection gain is at least its signed
projection gain because x_A(u)>=0. Thus the same distance lower bound
applies to all such signed planes.

If the unsigned partition is A/B but the signed plane differs from U_A,
some block contains two coordinates with opposite signs in its signed
line. That block has size at least2, even when the other block is a
singleton. The distance is at least sqrt(2) times the value of that
block, hence at least sqrt(2)u>=sqrt(2)/4. If both blocks have constant
signs, their plane is U_A itself and is not an exception.

Consequently every body pair plane U!=U_A satisfies

    dist(x,U)>=g_k:=1/(4sqrt(2))-sqrt(k)rho>0 on D_A.      (3)

The positivity gate is1/32>k/4096, true for all k<=5. The only body
lines in U_A are v and w: a body vector constant on each block has
block values+/-1, and its first coordinate fixes the A value+1.
Therefore no pair avoiding v, or avoiding w, spans U_A. All pair planes
used in the complementary-event estimates satisfy(3).

This says nothing about a uniform distance for x-a v or x-a w as a
varies. Such a claim is false, as043 demonstrated. The actual parameter
density must remain in the third-distance integral.

## 3. Complete positive branch bounds, with all Gram factors

Use the same actual partition as043:

    E={|a_v|>=1/2,|a_w|>=1/16},
    E_v={|a_v|<1/2},
    E_w={|a_v|>=1/2,|a_w|<1/16}.

Its densities P_s and B_s=B_v,s+B_w,s are nonnegative, sum to p_s,
and are each at most p_s. Enlarging E_w to |a_w|<1/16 is permitted for
the upper pointwise estimate only; the actual disjoint B_s is retained
when using its L1 bound.

The two selected columns have Gram determinant

    36-(k-m)^2=4km, J2=2sqrt(km).

By033, Q4/J2=R2,min(k,m)/(2pi)^4. Capping the two selected densities
and using genuine two-column coarea, exactly as in043, gives

    P_s<=R2,min(k,m)/(16pi^6) * s^-2/
                    [(1/4+s^2)(1/256+s^2)]
        <=D_big,k s^-2,
    D_big,k=64R2,min(k,m)/pi^6.                           (4)

The Gram root cancels only after its correct coarea factor is included.
This is why no k-dependent J2 remains in the last coefficient.

For a retained small parameter b=v or w, each partial quotient after
k'=0,1,2 non-b deletions has actual Gram32I-(P b)(P b)^T, between26I
and32I. The remaining parameter counts31,30,29 produce992s^2,960s^2,
928s^2. Each removed column has its real Gram-Schmidt coarea factor.
The final dimension3 quotient omits b; the proved symbol comparison
psi-|b.xi|>=(9/10)psi gives peak loss(10/9)^3. It is not replaced by
the original Q3 without that loss, nor by a rank4/Q2 quotient.

Keep the third denominator as a function of the retained parameter.
The exact scalar Cauchy average from043 gives, for original z!=0,

    integral_R f_s(a)/(|z-a b'|^2+928s^2) da
       <=(1+|b'|/sqrt(928))/|z|^2.

It includes all isolated distance zeros. Here |b'|<=sqrt(6), and(3)
bounds the original z distance. The full restricted ordered-triple
coefficient is C_rem=44743079/(2970pi^6), the same for any omitted body
line by symmetry. All26970 ordered triples are retained for each bound.

With F=1+sqrt(6/928), define

    L_k=(10/9)^3 C_rem F/g_k^2,
    M_v,k=L_k/(h_v0^2 h_v1^2),
    M_w,k=L_k/(h_w0^2 h_w1^2), M_k=M_v,k+M_w,k.           (5)

Then B_s<=M_k on D_A. All cancellation events and large parameters in
other directions remain paid. The arguments establishing pointwise
versions in043 apply unchanged: all quotient symbols have the stated
positive Gram lower bound, and the two-selected-column remaining full
frame has Gram32I-vv^T-ww^T positive definite (at least20I suffices for
all k). There is no assumption of a tight remaining frame.

## 4. The true joint marginal for each B size

Choose one fixed j in B. The original full marginal is always

    p2_s(1,z)=(1/2)f_(16s)((1+z)/2)f_(16s)((1-z)/2).

The remaining four-coordinate positive orthant has exactly1/16 of its
mass. This is a statement about p_s, not either branch. On the actual
zero-sum B residual, the j-coordinate functional has norm
sqrt((m-1)/m), including the value0 when m=1. Therefore

    L_m=1/4-sqrt((m-1)/m)rho,
    U_m=3/4+sqrt((m-1)/m)rho,
    integral_DA p_s dy<=C_m s^2,
    C_m=(128/pi^2)[F0(U_m)-F0(L_m)],
    F0(z)=z/[2(1-z^2)]+log((1+z)/(1-z))/4.               (6)

These endpoints lie strictly between0 and1. For m=1 the mean_B is the
coordinate itself, so its range is exactly[1/4,3/4] and no residual
margin is added. There is no multiplicity m in(6): choose a single j.

Since actual P_s,B_s<=p_s, (4)-(6), interpolation and Minkowski give

    ||S_epsilon||_L^(6/5)(D_A)<=U_k,
    U_k=(3/7)D_big,k^(1/6)C_m^(5/6)epsilon^(7/3)
         +(3/8)M_k^(1/6)C_m^(5/6)epsilon^(8/3).          (7)

The same proofs cover k=1 and k=5; no empty-block quotient or density
factor is introduced. Only the zero-dimensional residual of a singleton
block is interpreted as0.

## 5. k3 source reuse and the31-domain union

For k=m=3, all formulas(2)-(7) reduce symbolically to043 and049: eta is
sqrt(3)rho, the squared unshifted line gaps are3/32 and3/128, J2=6,
and the B residual norm squared is2/3. The actual k3 scalar coefficients
and U3 are read from complete saved intervals. They are not recomputed.
Only k=1,2,4,5 have new geometric/scalar configurations.

Permuting the five free physical coordinates maps any A of size k
containing1 to the representative A={1,...,k}, preserving the full
original frame, positive face measure and D_A. The exact counts are

    n_k=binom(5,k-1)=(1,5,10,10,5), sum n_k=31.

For the common nonnegative function S_epsilon and q=6/5,

    integral_(union_A D_A) S_epsilon^q
       <=sum_A integral_DA S_epsilon^q
       <=sum_(k=1)^5 n_k U_k^q,

    ||S_epsilon||_Lq(union_A D_A)
       <=[sum_(k=1)^5 n_k U_k^(6/5)]^(5/6).              (8)

Overlaps are permitted and counted positively. The31 domains are not
asserted disjoint, and no k versus6-k identification removes weights:
the first coordinate is in A, and its geometric role is not symmetric
with B. Formula(8) covers exactly the31 specified middle domains. It
does not cover the spatial endpoints, longer times, the whole face,
the global g condition or the full I/g cost.
