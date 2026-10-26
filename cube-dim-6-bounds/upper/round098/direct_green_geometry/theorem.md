# 098: an exact time-first marginal identity and one paid local Green region

All laws in this note are the original equal-weight body32 Cauchy model:
V={(1,s2,...,s6):si=+/-1}, psi(xi)=sum_v |v.xi|,
p_t=(2pi)^(-6) Fourier inverse of exp(-t psi), and G=int_0^infinity p_t dt.
The scalar Cauchy density is c_t(a)=t/[pi(a^2+t^2)]. All parameters have the
same time t. All projections and quotient densities use orthonormal Lebesgue
measure. Face integration uses dy on x=(1,y), y in [0,1]^5.

This is an author proof submission, conditional on the read-only accepted
033/037/043 source interfaces. The four arithmetic configurations do not
replace review of the proof. The model is unchanged; no global g or complete
competitive I is claimed.

## 1. A true four-dimensional time integral equals a three-dimensional marginal

Fix two independent body directions v,w, U=span(v,w), W=U^perp. Their
projections vanish. The remaining 30 projected columns b_j=P_W v_j are
nonzero, span W, and satisfy sum_j b_j b_j^T=32 I_W. Denote their actual
four-dimensional common-time density by q_t, with symbol psi restricted
to W. For a unit e in W, let Q3(e) be the time-one density at zero of the
orthogonal projection of this law onto E_e=W intersect e^perp. Then, for
**every** z!=0 in W,

    H(z):=int_0^infinity t^2 q_t(z)dt = Q3(z/|z|)/(2|z|).       (1)

This does not assert a pointwise isotropic Cauchy bound for q_t.

Here are the pointwise and finiteness details. The positive deletion proof,
using |sum b_j a_j|^2<=32 sum a_j^2, gives at t=1

    q_1(x) <= [32/(pi(|x|^2+960))]
                 sum_j |b_j|^(-1) q3_j,1(P_(b_j^perp)x)
           <= C/(|x|^2+960),                                (2)

where C=(32/pi)sum_j |b_j|^(-1)q3_j,1(0)<infinity. Each q3_j is a
genuine three-dimensional projected law: the projected remaining columns
span that three-space. Coercivity gives its finite continuous Fourier
density and the peak bound q3_j,1(y)<=q3_j,1(0). The original q_1 is
continuous as well. The test-function deletion inequality therefore holds
pointwise, as in the accepted positive coarea interface.

The marginal identity

    int_R q_1(se+y)ds = q3_e,1(y)

first holds almost everywhere in y in E_e by Tonelli. Its left side is
continuous in y: since y perpendicular e, (2) is dominated by
C/(s^2+960), uniformly for y, and q_1 is continuous. The right side is
continuous by its coercive Fourier integral. Equality thus holds at every
y, in particular y=0. Evenness now gives

    int_0^infinity q_1(se)ds = Q3(e)/2.

Finally q_t(re)=t^(-4)q_1((r/t)e). The positive change s=r/t yields

    int_0^infinity t^2 q_t(re)dt
        = (1/r)int_0^infinity q_1(se)ds.

This proves (1) for every direction, including a projected spectral
direction at which q_t(re) can have a short-time spike. No Fourier Abel
limit, interchange of oscillatory integrals, or equality merely almost
everywhere is needed. The bound (2) is used only to justify a continuous
line trace, not as the new numerical coefficient.

Three-dimensional Fourier inversion and polar integration also give

    Q3(e) = (1/(4pi^3)) int_(S(W) intersect e^perp) psi(theta)^(-3) dS2,
    H(z)  = (1/(8pi^3|z|)) int_(S(W) intersect z^perp) psi(theta)^(-3) dS2.
                                                               (3)

Thus the suggested positive spherical expression has the stated constant.
It follows here from an ordinary positive marginal identity. At z=0 the
time integral diverges; no finite value is claimed there.

## 2. A uniform coefficient from exact even moments

For independent uniform signs eps_i put S=sum eps_i xi_i. The 32 unoriented
body lines give psi(xi)=32 E|S|. The exact moments are

    E S^2=|xi|^2,
    E S^4=3|xi|^4-2sum xi_i^4 <=3|xi|^4.

Holder, applied to |S|^(2/3)|S|^(4/3), gives
(E S^2)^(3/2)<=E|S| (E S^4)^(1/2). Consequently

    psi(xi)>=c|xi|, c=32/sqrt(3).                            (4)

The zero case is immediate. Restriction to every actual subspace preserves
(4). It follows from the three-dimensional radial Fourier integral that

    Q3(e)<=1/(pi^2 c^3),
    H(z)<=1/(2pi^2 c^3 |z|).                                (5)

This is an integrated estimate. It does not imply a pointwise bound on q_t.

One optional directional version retains a saved exact peak. For
v=(1,1,1,1,1,1), w=(1,1,1,-1,-1,-1), and
a=(1,1,-1,1,1,-1), put e0=P_W a/|P_W a|. The exact Gram determinants are
J2^2=36 and J3^2=192. Hence 033 gives

    Q3(e0)=sqrt(192)*(707/6614784)/(2pi)^3.

For ||e-e0||<=1/64 take the plane rotation R with R e0=e. Since
psi(h)<=32|h|, ||R-I||=||e-e0||, and the accepted elementary bound
psi(xi)>=10|xi| holds, psi(R xi)>=(19/20)psi(xi). An orthonormal change
of variables in the three-dimensional peak integral proves

    Q3(e)<=(20/19)^3 Q3(e0).                                (6)

No new fan or projected peak was evaluated. This optional cone is not
needed for the full-direction local region below.

## 3. Keep the true two-large-parameter source before integrating time

Use precisely v,w above and the actual exhaustive event split from 043:

    E={|a_v|>=1/2, |a_w|>=1/16},
    E_v={|a_v|<1/2},
    E_w={|a_v|>=1/2, |a_w|<1/16}.

Let P_t be the density on E and B_t the density on its actual disjoint
complement. Thus p_t=P_t+B_t. On E the two selected Cauchy factors are at
most t^2/[pi^2(1/4+t^2)(1/256+t^2)]. Push this positive bound through
the full source fiber. The free two-parameter coarea is 1/J2=1/6, and
the other 30 columns give the actual q_t. Therefore, for all x and t>0,

    P_t(x) <= t^2 q_t(P_W x)
                 /[6pi^2(1/4+t^2)(1/256+t^2)].              (7)

In particular the tangential denominators can be retained in the positive
time integral. If z=r e, the same change of variable as in (1) gives the
more informative expression

 int_0^infinity P_t(x)dt
 <= (1/(6pi^2 r)) int_0^infinity
       q_1(s e)/[(1/4+r^2/s^2)(1/256+r^2/s^2)] ds.          (8)

No tangential/normal independence has been assumed. The true tangential
parameter values depend on the whole source fiber; (7) used valid event
caps before freeing these two parameters. For the present finite estimate
drop only the two positive r^2/s^2 terms. Equations (1) and (5) give

    int_0^infinity P_t(x)dt <= K/r,
    K=1024/(12pi^4 c^3)=sqrt(3)/(128pi^4)<1/7000.           (9)

This covers both short and large positive times. It is stronger than a
per-time Q4 peak near r=0, and it allows the true anisotropic short-time
spikes of q_t.

## 4. All complementary source events are paid on one explicit outer tube

Let D be the 043 domain: x=(1,y), y in [0,1]^5,
rho=dist(x,U)<1/64 and u=mean_(4,5,6)x in [1/4,3/4]. Set

    eta=sqrt(3)/64,
    hv0=sqrt(3)/2-eta,       hv1=sqrt(3/32)-eta,
    hw0=15sqrt(3)/16-eta,    hw1=sqrt(3/128)-eta,
    g=1/(4sqrt(2))-eta,
    Crem=44743079/(2970pi^6),
    Krem=(1000/729)Crem(1+sqrt(6/928)).

All five gaps are positive. The accepted actual-fiber bound 043, retaining
the time terms as in 074, is

 B_t <= Krem/(g^2+928t^2)
       [1/((hv0^2+992t^2)(hv1^2+960t^2))
        +1/((hw0^2+992t^2)(hw1^2+960t^2))].                 (10)

It includes all other large parameters, cancellation events, missing-column
quotient loss, and isolated zeros of a shifted third distance. No branch
reflection symmetry is used. Dropping the last four positive time terms
and integrating only the retained resolvent gives the new coarse all-time
payment

 int_0^infinity B_t dt <= B0
   =pi Krem/(2sqrt(928)g)
          [1/(hv0^2 hv1^2)+1/(hw0^2 hw1^2)] <512.           (11)

The actual Arb enclosure has B0 about343.378053668633. This is one new
analytic integral, not a replay/refinement of 074's time rectangles or
092/096's excluded information relaxations. Equations (9)-(11) prove

    G(x)<=1/(7000rho)+512 on D, rho>0.                    (12)

Every positive time and every source event is included. At rho=0 the
right side may be defined as infinity.

## 5. A strictly better upper expression on a finite new local region

Freeze delta=2^(-24) and define

    Q={x in D: 31/64<=u<=33/64, 0<rho<=delta}.

This is a nonempty positive-measure region. With x(u)=(1,1,1,u,u,u),
|x-x(u)|<=sqrt(3)rho. Even replacing rho by the larger1/64 gives, on Q,

    |x|^2 <= [sqrt(3(1+(33/64)^2))+eta]^2 <4,
    dist(x,span v)^2
        <= [sqrt((3/2)(1-31/64)^2)+eta]^2 <1/2.             (13)

The saved 037 pair-resolved all-time envelope Phi contains its positive
ordered (v,w) term M3/[|x|^2 dist(x,span v)^2 rho]. No term is subtracted.
Its frozen dyadic lower endpoint proves M3/2>1/5000. Thus

    Phi(x)>=1/(5000rho),
    G(x)<=1/(7000rho)+512
        <=(24855/28672)Phi(x), x in Q.                     (14)

The exact ratio is5000(1/7000+512/2^24)<1. This is a strict improvement
of a previous upper **expression** on Q. It does not assert that the true
G decreased, or that the original model cost has been excluded. The very
small delta is fixed once; no radius scan followed the result.

## 6. A full-time local spatial q budget

Only now take the spatial L^(6/5) norm. Write x_A=alpha1_A+e,
x_B=u1_B+z, sum e=sum z=0. The face x1=1 sets alpha=1-e1, and
rho^2=|e|^2+|z|^2. The face Jacobian from orthonormal (e,z,u) coordinates
is3: the A map (e2-e1,e3-e1) has Jacobian sqrt(3), and the B mean map
has Jacobian sqrt(3). The constraints x2,x3<=1 imply e1>=e2,e3, a
maximal-coordinate sector with exactly1/3 of radial measure in the A
zero-sum plane. The other face constraints only shrink the domain.

Let v4=pi^2/2 and q=6/5. The u interval has length1/32, so

    |Q| <= (v4/32)delta^4,
    int_Q rho^(-q)dy <= (5v4/112)delta^(14/5).              (15)

Spatial Minkowski applied to the already integrated (12) proves

    ||G||_(Lq(Q)) <= B_Q
      =(1/7000)[(5v4/112)delta^(14/5)]^(5/6)
          +512[(v4/32)delta^4]^(5/6),
    int_Q G^q dy <= B_Q^(6/5) =: I_Q,U <10^(-23).          (16)

The actual outward enclosure is B_Q about6.511925611508e-22 and
I_Q,U about3.770954714088e-26. These are deterministic upper budgets,
not sampled values or estimates of the actual local integral. Q lies in
the true middle tube; the other pair planes are separated by the 043 gap
g>1/64, so this particular small region meets no other radius1/64 tube.
No assumption about disjointness of the full31-tube family is made.

## 7. Remaining operator and cost gaps

For the same body32 model the objective is

    S=64(g_L/6)^(6/5)-I_U,
    I_U >= int_[0,1]^5 G(1,y)^(6/5)dy.

The candidate g0=1/65000 is still not a global certificate. Conditional
on it, the whole-face norm gate is B<4/48750. This note pays the full G
only on Q and supplies a new general normal-potential identity. To insert
(16) into a complete ledger, give Q its own spatial region and certify the
complementary face integral without omission or duplicated allocations.
No complete competitive complement or matching global g_L is supplied.

The framework still targets Lebesgue measure, centered axis-parallel cubes,
one common radius in all coordinates, general L1 input, and every positive
radius. Nothing here restricts that target to atoms or a finite radius set;
the missing global certificates remain essential. The established C6 upper
64 and accepted strict lower2073015624830/10^12 are unchanged.
