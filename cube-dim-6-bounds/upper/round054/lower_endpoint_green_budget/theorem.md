# Full positive-time Green budget on every unsigned lower endpoint tube

Use only the original body32 equal-scale Cauchy frame
V={(1,s2,...,s6):si=+/-1}. The face measure is dy on
F={(1,y):y in[0,1]^5}, all distances are Euclidean in R6, and
G(x)=integral_0^infinity p_t(x)dt. Fix rho=1/64 and q=6/5.
For A containing1, 1<=k=|A|<=5, let B=A complement, m=6-k,
U_A=span(1_A,1_B), and

    E_A={x in F: dist(x,U_A)<rho, u=mean_B(x) in[0,1/4]}.

This is the full-time lower endpoint tube, not a short-time kernel budget.
The envelope input is C037 theorem(7), independently accepted by F037:

    G(x)<=Phi(x)=sum_(ordered v!=w) M_class(v,w)/
                 [|x|^2 dist(x,span v)^2 dist(x,span(v,w))].       (1)

It holds off S2 and with +infinity at singular denominators. All 992
ordered terms occur, with class counts (192,480,320). The positive M_j
are read as complete saved intervals from C037; their defining integrals
are not reevaluated in this task.

## 1. Uniform norm and line bounds on every E_A

Write x_A=alpha1_A+e, x_B=u1_B+z, with both residuals zero-sum.
As x1=1, e1=1-alpha. For x_A(u)=1_A+u1_B,

    h=x-x_A(u)=(e-e1 1_A,z),
    |h|^2=dist(x,U_A)^2+k e1^2<=k dist(x,U_A)^2,

because e1^2<=(k-1)|e|^2/k. The same formula includes k1, when e=0.
Set eta_k=sqrt(k)rho. Then

    |x|>=r0,k=sqrt(k)(1-rho)>0.                            (2)

For x_A(u), nonnegative coordinates imply that every signed body line
has projection gain at most that of ones. Its exact minimum distance is
sqrt(km/6)(1-u). Since 0<=u<=1/4, and distance to any subspace is
1-Lipschitz, every body line has distance at least

    d0,k=(3/4)sqrt(km/6)-sqrt(k)rho>0.                    (3)

Positivity follows exactly from (3m/32)>rho^2 for m>=1; no sampling is
used. The zero B coordinates at u0 may produce tied nearest lines, but
not a smaller distance. Thus singleton blocks and partial zero
coordinates cause no exception.

## 2. All incident pair planes, without full direction enumeration

Any body pair v,w has agreement block C containing1 and nonempty
complement D. Its plane has orthogonal support directions v1_C and
v1_D. Suppose that it contains the endpoint 1_A. A nonzero coefficient
on either direction produces a nonzero coordinate at every index in
its support, since body signs never vanish. Because A is nonempty and
proper, exactly one coefficient can be nonzero. Therefore its support
must equal A. Since 1 is in A and C, it must be C=A; the A signs are
constant, while the B signs are unrestricted. Conversely every plane

    W_sigma=span(1_A,sigma_B), sigma_B in{+/-1}^B,

contains1_A and is a body pair plane, represented by
v=(1_A,sigma_B) and w=(1_A,-sigma_B).

Two choices sigma and tau give the same B line exactly when tau=+/-sigma.
Hence there are exactly2^(m-1) incident planes. Each contains precisely
two body lines: a body vector in the plane has A value+1 (its first
coordinate is1), and B coefficients must be either+sigma or-sigma.
Thus each plane contributes exactly two ordered pairs, for a total
2^m. Their Hamming separation is m, hence all belong to
j=min(k,m). This proves completeness, including k1 and m1, without
enumerating 496 planes or 992 ordered pairs.

For x>=0, the incident signed plane has

    dist(x,W_sigma)^2
      =|x|^2-(sum_A x_i)^2/k-(sum_B sigma_i x_i)^2/m
      >=|x|^2-(sum_A x_i)^2/k-(sum_B x_i)^2/m
      =dist(x,U_A)^2.                                    (4)

Mixed B signs are retained in (4), not discarded or assumed irrelevant.

## 3. Every nonincident signed plane stays uniformly separated

At x_A(u), if its unsigned partition differs from A/B, one of its blocks
contains both an A coordinate of value1 and a B coordinate of valueu.
The squared error of those two coordinates relative to any common
constant is at least(1-u)^2/2. Signed projection gain at a nonnegative
point is no larger than unsigned projection gain on the same partition.
Consequently all these signed planes have distance at least
(1-u)/sqrt(2)>=3/(4sqrt(2)).

If the partition is A/B, a nonincident plane must have nonconstant signs
on A: all constant-A sign choices with arbitrary B signs were included
in section2. For opposite A signs, two A coordinates of value1 have
squared residual at least2 relative to a signed constant. Thus distance
is at least sqrt(2), which is stronger than3/(4sqrt(2)). When k1 such a
sign pattern is impossible, so the case is empty. These alternatives
cover all nonincident signed planes, including those mixed on both
blocks. Subtracting |h| gives

    dist(x,W)>=g_k=3/(4sqrt(2))-sqrt(k)rho>0.              (5)

The exact positivity gate is9/32>k rho^2 for k<=5.

## 4. Full-time pointwise envelope, with all background terms charged

Let N=(192,480,320), j=min(k,m), c_k=2^m, and set

    I_k=c_k M_j/(r0,k^2 d0,k^2),
    B_k=[sum_(ell=1)^3 (N_ell-c_k 1_{ell=j}) M_ell]
                                      /(r0,k^2 d0,k^2 g_k).     (6)

Every coefficient N_ell-c_k 1_{ell=j} is a nonnegative integer, and the
background counts sum to992-c_k. Combining (1)-(5), with
s=dist(x,U_A), gives the complete bound

    G(x)<=Phi(x)<=I_k/s+B_k on E_A.                       (7)

The numerical implementation uses the nonnegative sum in (6); it does
not subtract two correlated interval evaluations of Mtotal. Incident
terms that are actually nonsingular at a given positive point are still
charged by I_k/s. This intentional enlargement is harmless.

All S2 intersections inside E_A are among the incident planes by (5).
The inequality may take the value+infinity on their null sets, as may
G. It remains a valid a.e. upper bound for the Lq integral. A finite
union of their face intersections has codimension4 and measure zero.
No actual asymptotic coefficient or lower bound for G is asserted.

## 5. True face measure, radial singular integral and boundary clipping

Use the F034 coordinates. The A residual e lies in the zero-sum space of
dimension k-1 and the face coordinate condition xi<=1 is exactly the
cone C_A={e1>=e_i}. The map from A deficits to e has Gram
I_(k-1)-11^T/k, determinant1/k, so its inverse face Jacobian is sqrt(k).
For B the mean u is not a unit coordinate: dx_B=sqrt(m)du dz.
Together dy=sqrt(km)du de dz. The distance s is exactly|(e,z)|.

The normal space has dimension(k-1)+(m-1)=4. Permuting A coordinates
preserves its radial measure and partitions it into the k maximal-
coordinate cones, up to zero-measure ties. Hence the cone fraction is
1/k. This includes k1 with zero-dimensional A residual and fraction1;
likewise m1 simply has no B residual. No orthant factor2^(-k+1) is used.

Keep 0<=u<=1/4 and the cone, and enlarge all other box restrictions.
In particular u0 may force additional B-coordinate clipping; deleting
it can only increase these nonnegative integrals. With v4=pi^2/2,

    V_k:=|E_A| upper
       =(1/4)sqrt(m/k)(pi^2/2)rho^4,
    J_k:=integral_EA s^(-q) upper
       =(1/4)sqrt(m/k)[2pi^2/(4-q)]rho^(4-q).             (8)

These are upper bounds, not claims of actual volume equalities. The
radial integral converges since q<4; endpoint clipping is fully allowed.
Minkowski and (7)-(8) give a full positive-time budget

    ||G||_Lq(E_A)<=U_k:=I_k J_k^(5/6)+B_k V_k^(5/6).      (9)

## 6. All31 domains, with overlaps and limitations

Permutations of the five free coordinates preserve the frame, G and
face measure, so there are n_k=binom(5,k-1) domains of type k. The exact
weights (1,5,10,10,5) sum to31. Positivity gives

    ||G||_Lq(union_A E_A)
       <=[sum_(k=1)^5 n_k U_k^(6/5)]^(5/6).              (10)

The domains may overlap; (10) does not assume disjointness or charge a
signed plane as a new physical domain. All 992 ordered terms were
included for each E_A before this positive union estimate.

This covers the low-mean tubes only, over all t>0. It does not pay the
high-mean endpoint u>3/4, the middle-domain long-time tail, outside-tube
face region, or any global g requirement. If the fixed upper formula
exceeds an intended target, that closes this formula as a target budget,
not the possibility of a sharper bound or the actual Green norm. No
radius scan, new kernel or parameter repair is part of this batch.
