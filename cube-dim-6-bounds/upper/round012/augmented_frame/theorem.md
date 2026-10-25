# Augmenting body32 by six coordinate axes

## New model and exact short-potential interface

Let V={(1,s2,...,s6):s_i=+/-1}. For beta>0 let
A_beta consist of these 32 columns and beta e1,...,beta e6, and put

    psi_beta(xi)=sum_(v in V)|v.xi|+beta sum_i |xi_i|.

If the 38 scalar variables are independent Cauchy variables of scale t,
their image by A_beta has characteristic function exp(-t psi_beta).
Thus its density p_beta,t is the convolution of the body32 density with
the six independent coordinate Cauchy densities of scale beta t. These
are symmetric conservative probability kernels. Convolution defines a
positive, mass-preserving, constant-preserving strongly continuous L1
contraction semigroup P_beta,t. Spanning, and hence a smooth fixed-time
density, follows already from the six added coordinate columns. In
particular psi_beta>=beta ||xi||_1 ensures Fourier integrability. No
negative weights or subtraction of stochastic kernels defines this model.

The generator is L_beta=-psi_beta(D). For T>0, Tonelli gives

    K_beta,T=integral_0^T p_beta,t dt >=0,
    ||K_beta,T||_1=T.

For compactly supported smooth phi, differentiation of the semigroup,
or its Fourier multiplier, gives distributionally

    L_beta K_beta,T=p_beta,T-delta_0.                  (1)

The order-one symbol gives p_beta,t(x)=t^-6 p_beta,1(x/t). Consequently
the usual spatial dilation obeys

    (K_beta,T)_r(x)=r^-6 K_beta,T(x/r)
                  =(1/r) integral_0^(rT) p_beta,s(x)ds,
    ((K_beta,T)_r*f)/T=(1/(rT)) integral_0^(rT)P_beta,s f ds. (2)

These are all positive-time averages of one fixed Markov semigroup.
Their weak (1,1) constant is one, for all nonnegative L1 functions. For
completeness, fix h>0, P=P_beta,h, f>=0 and lambda>0, and define
R0=0, Rn=(f-lambda+P R_(n-1))_+, Dn={Rn>0}. Positivity implies
Rn>=R_(n-1) and Rn<=f+P R_(n-1) in L1. The inclusion
Dn subset {f+P R_(n-1)>lambda} and Chebyshev show that Dn has finite
measure. Integrating on Dn,

    integral_Dn(f-lambda)
       =integral Rn-integral_Dn P R_(n-1)
       >=integral Rn-integral R_(n-1)>=0.

Also Rn>=sum_(j=0)^(n-1)P^j f-n lambda. Therefore the union of all
discrete-average bad sets is contained in the increasing union of Dn,
of measure at most ||f||_1/lambda. For finitely many positive rational
times take mesh h dividing them. Strong L1 continuity bounds the error
of each left Riemann average by sup_(0<=s<=h)||P_beta,s f-f||_1,
which tends to zero. The finite maxima converge in L1 and retain the
weak bound by Chebyshev and a vanishing threshold loss. Exhaust the
rational times and use almost-everywhere continuity of the indefinite
time integral (Fubini applies since integral_0^R ||P_beta,s f||_1 ds
is finite). This proves the all-positive-time bound with constant one.

If K_beta,T>=m>0 almost everywhere on [-1,1]^6, then for every radius r
the centered cube average of |f| is bounded by
((K_beta,T)_r*|f|)/(64m). Rational radii and the usual continuity of
cube averages give the common all-radius exceptional-set formulation.
Equation (2) and the proved maximal bound therefore imply, for general
L1 input and Lebesgue measure,

    C6 <= T/(64m).                                     (3)

This is a newly checked interface for the augmented semigroup. It does
not inherit an old body32 value of g, I, m, or an old singular-set formula.
For the registered T=2/45, the necessary threshold for improvement below
64 is m>T/4096=1/92160. The sufficient target
m0=T/4000=1/90000 would give 125/2 only if certified globally.

## Why the old exact fan geometry survives

The original 21 exact simplices partition the Fourier chamber
xi1>=xi2>=...>=xi6>=0 up to boundaries, after taking xi1=1 as section.
Their 46080 signed-permutation images cover all R6 up to null boundaries.
This geometry is the existing rational certificate and its independent
exact reconstruction, not a new triangulation inferred from decimals.
Every body linear form has constant weak sign on each simplex. Each
coordinate form is nonnegative there, so sum_i |xi_i| is the linear
function sum_i xi_i throughout that whole cone. Therefore on a cone
with integer rays r_j,

    q_beta,j=psi_beta(r_j)=q_body,j+beta sum_i r_(j,i).

The rays and determinant are unchanged; the new alpha is
|det(r1,...,r6)|/product_j q_beta,j. Signed-permutation invariance of
psi_beta makes the same full group valid. The finite audit verifies every
ray occurrence against the original rational certificate and the explicit
38-column symbol, not only the displayed rate-update formula.

For b_j=x.r_j/q_beta,j, Fourier integration on a cone thus gives exactly
the same universal H formula with the new q and alpha:

    H_beta,T=(2pi)^-6 sum_R alpha_R DD5(atan(./T);b_R).

This identity depends on linearity on cones, not on the number of frame
columns. Confluent divided differences use derivatives divided by the
factorial, including zero nodes. As T decreases to zero, the full
Laurent sum retains its constant coefficient and all possible
T^-1,T^-3,T^-5 principal coefficients. They must cancel exactly at any
finite Green point; no individual cone term can be discarded.

## General finite-rank regularity, including dependent triples

Let S2_beta be the union of spans of at most two of the 38 columns.
Distinct triples need not be independent: if body vectors differ in
one coordinate, their difference is twice that coordinate axis.
Nothing below assumes otherwise.

On a compact set E away from S2_beta choose eta=dist(E,S2_beta)>0 and
delta=eta/[2 sum_j ||a_j||]. If x=sum a_j u_j is in E, the columns
whose parameters obey |u_j|>=delta have rank at least three; otherwise
the sum of the remaining terms places x within eta/2 of a span of at
most two columns. Select any independent triple of large parameters
and extend it to a six-column basis B. Solving those six variables gives
the exact positive density formula with prefactor 1/|det B|. Its three
large Cauchy factors are at most t/(pi delta^2), and its three other
basis factors at most 1/(pi t). Integrating the remaining 32 probability
densities and summing over finitely many independent triples gives

    p_beta,t(x)<=C_E,   x in E, 0<t<=1.                (4)

This is a rank-selection argument: dependent triples are never selected.
To make the density inequality formal, multiply the parameter-space
event cover by a nonnegative test function phi(A_beta u), integrate the
38 independent Cauchy densities, and on each selected event change the
six basis variables to x. Its Jacobian is |det B|. The event-restricted
pushforward density is then bounded by the displayed three-large-factor
estimate, while the remaining 32 scalar densities integrate to one.
First use the eta/4 open neighborhood of E for these test functions; it is
still at distance >=3eta/4 while the discarded remainder is <eta/2,
so the event cover still holds. The resulting a.e. density inequality
on that open neighborhood extends to E by fixed-time continuity. No
conditioning on a zero-probability fiber has been used.

For large t, Fourier integrability and scaling give p_beta,t<=C t^-6.
Dominated convergence therefore proves G_beta=integral_0^infinity
p_beta,t finite and continuous off S2_beta, and there
K_beta,T=G_beta-H_beta,T. This subtraction is used only at such points.
For the two registered strict positive points, six distinct nonzero
absolute coordinates exclude all three possible types of two-column
span: body/body has at most two absolute values; body/axis has equal
absolute values in the five unmodified coordinates; axis/axis has at
most two nonzero coordinates. This argument is only a membership check,
not a closest-plane geometry enumeration.

Thus a strictly subthreshold certified K value at either face point
extends to a relative neighborhood in the cube, with positive volume,
and really obstructs the essential lower bound. Survival at one or two
points provides only necessary conditions. It supplies no global lower
certificate and no C6 upper bound.
