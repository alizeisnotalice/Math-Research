# F028: exact all-time acceptance and unrestricted probability mixtures

This file states the proof being audited. Its numerical premises must be
verified by the independent root checker and the two native H payments;
this file alone is not a completed all-time certificate.

Fix body32 and x0=(1,1,1,9/10,27/50,0). Set

    c=1/4100, ell=3141592653/10^9, u=3141592654/10^9,
    K(T)=integral_0^T p_t(x0)dt,
    alpha=(pi/ell)^6,
    F(T)=cT-alpha K(T).

The point is off S2 (its four coordinate absolute levels cannot occur in
a body-pair plane), so the independently proved finite-frame regularity
gives G(x0)<infinity and K(T) finite. F023 independently evaluated this G
as Gamma/(128*pi^5), retaining all zero-node Laurent contributions and
checking their exact cancellation.

## Rational identity without a new global numerator

The B023 density table represents R(z)=sum d_(s,j)/(z+s)^j with positive
distinct s. Its true denominator is the monic
D(z)=product (z+s)^m_s, of degree417; its 359 highest local coefficients
are nonzero. B024 stores a proper N of degree414.

For each s, expand E_s(z)=D(z)/(z+s)^m locally at z=-s by multiplying the
jets of its other linear factors. Then expand N/E_s to order m-1 by
ordinary power-series division. Matching the resulting 417 principal
coefficients with B023 proves N/D=R: their difference is a rational
function without finite poles and vanishes at infinity, hence is zero.
The nonzero highest coefficient proves N(-s)!=0 and gcd(N,D)=1. No
global denominator expansion or new PF numerator synthesis is required.

Let Z be the saved primitive integer polynomial, a>0 its clearing scale,
and k=64ell^6 c. To verify P=N-kD=aZ independently, check that N-aZ
vanishes to multiplicity m_s at every -s. These are 417 zeros counted
with multiplicity. Since deg(N-aZ)=417, comparison of its leading
coefficient with k determines N-aZ=kD exactly. This establishes both
the function identity and the positive sign of the clearing operation.

D>0 for z>0 because every factor z+s is positive. Thus

    p_T(x0)=N(T^2)/(64pi^6 D(T^2)),
    F'(T)=c-alpha p_T(x0)=-P(T^2)/(64ell^6 D(T^2)).

The negative sign in this derivative is essential. F(0+)=0 because
K(T)->0, and F(T)->+infinity as T->infinity because K(T)<=G(x0)<infinity.
The alternative potential alpha K-cT has all three signs reversed;
mixing the two conventions invalidates the extrema argument.

## Complete positive critical roots

Verify the saved exact squarefree factorization reconstructs Z, each
factor's gcd with its derivative is one, and the squarefree part S has
nonzero constant term. Verify its coefficient sign-variation count is2.
Two disjoint positive rational intervals with opposite nonzero endpoint
signs supply at least two distinct positive roots. Descartes' rule permits
at most two counted with multiplicity. Hence the two intervals contain
all positive roots, exactly one each, and both are simple. The explicit
factor ownership and squarefree checks also verify their stated
multiplicities; no Sturm or new root isolation is needed.

The map T->T^2 is strictly increasing with nonzero derivative for T>0.
Exact inequalities a_T^2<=a_z<b_z<=b_T^2 transfer the two complete z
brackets to the locked positive T brackets and preserve multiplicities.
There are no other positive critical points of F, because its denominator
never vanishes there. The initial sign P(0)<0 and its two simple positive
roots imply signs negative, positive, negative on the complementary
positive intervals; F' has the opposite signs.

## Paying the full critical brackets independently

For a locked bracket [a,b] containing a critical point, K is nonnegative
and increasing. If a certified K_upper(b) satisfies

    c*a-(u/ell)^6*K_upper(b)>0,                         (1)

then for every T in that whole bracket,

    F(T)=cT-alpha K(T)
        >=c*a-(u/ell)^6*K_upper(b)>0.

The independent K enclosure is constructed as G-H_b. The G coefficient
is F023's already audited full rational Gamma, with an independently
computed rigorous pi enclosure. H_b is evaluated by the existing C010
native512 confluent-DD engine with all zero/repeated nodes and all full
group weights. This needs no new PF evaluation and no full G recomputation.
B028 PF bounds are read only for a crosscheck, not used in (1).

Once both complete critical brackets are paid, F has fixed derivative
sign on each remaining interval. On each bounded gap its minimum is at an
endpoint. On the first gap the limiting endpoint value is F(0)=0; on the
last, the infinity limit is +infinity. The bracket payments and these
limits therefore prove F(T)>=0 for every T>0. In fact the initial positive
derivative and strict payments give strict positivity away from zero.
Consequently K(T)/T<=c/alpha<c for all positive T. The conservative claimed
uniform bound remains c=1/4100.

These are roots of the rational comparison polynomial P using ell, not
necessarily the roots of the true equation p_T=c. No assumption of
unimodality of K(T)/T or p_T enters the argument.

## Arbitrary probability time mixtures: no inverse-time moment

Let mu be any probability measure on (0,infinity), with no moment
restriction, and define

    Q_mu(y)=integral K_T(y)/T dmu(T).

Each integrand is a nonnegative kernel of mass1, so Tonelli gives
integral Q_mu=1. For spatial mass-preserving dilation,

    r^-6 K_T(y/r)/T=K_(rT)(y)/(rT).

Thus every dilation convolution of Q_mu with |f| is an average of the
same Markov semigroup's prefix-time averages, and is bounded by its
ergodic maximal function. The previously independently proved positive
Markov weak-(1,1) theorem applies on infinite Lebesgue space to general
L1 and all positive scales. Neither swapping a supremum into convolution
nor separately optimizing each coordinate is involved.

The all-time point bound gives Q_mu(x0)<=c. To make this an essential
minimum obstruction, take a compact neighborhood E of x0 disjoint from
S2. The positive-parameter rank argument gives a common p_t(y)<=C_E for
0<t<=1 and y in E. A fixed six-column basis gives p_t(y)<=C t^-6 for
t>=1. Therefore G(y)<=G_E<infinity uniformly on E and

    K_T(y)/T<=C_E   (T<=1),
    K_T(y)/T<=G_E   (T>=1).

For each fixed T>0 the kernel is continuous on E. The displayed bound
is integrable against every probability mu, even if its support is
unbounded or accumulates at zero. Dominated convergence makes Q_mu
continuous near x0. Every neighborhood of this cube-boundary point
intersects the cube interior in positive six-dimensional measure, so

    essinf_cube Q_mu <= Q_mu(x0) <=1/4100.

Hence the standard domination model cost 1/(64 essinf Q_mu) is at least
1025/16. This concerns all normalized positive time mixtures of this
fixed body32 semigroup, not mixed frames or the actual value of C6.

No assumption integral T^-1 dmu<infinity is needed anywhere above.
If one separately writes LQ_mu as a positive measure minus a scalar
multiple of delta0, that finite-source decomposition would require
integral T^-1 dmu<infinity. Without it one can use the compensated test
function expression, but must not split two infinite terms. This audit
uses no such source decomposition. Likewise mixing raw K_T against nu
and normalizing mass would require 0<integral T dnu<infinity; that is a
different parameterization from the probability mixture above.

This argument has no implication for a different Cauchy frame, signed
time combinations, or other upper-bound mechanisms. The C6 endpoint is
unchanged by a successful model-family obstruction.
