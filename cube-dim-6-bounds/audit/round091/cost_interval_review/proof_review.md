# A091 independent core proof review

The fixed model is body32_unit_cauchy_v1 in dimension6, with 32 independent
unit-weight Cauchy(t) parameters at the same t. We use q=6/5,
F=[0,1]^5, G(x)=integral_0^infinity p_t(x)dt,
g=essinf_[-1,1]^6 G and I=integral_F G(1,y)^qdy. The notation U is a
constant-truncation model cost, while S denotes a surplus. No n5, different
alpha, different frame, or time-mixture endpoint is substituted.

The parent-supplied 091 manifest and every registered input hash are pinned.
The new guidance was read in full and its exact hash is recorded. This
review accepts previous 036/057/052 and the final 061/065 certificates as
explicit dependencies, without rerunning their peaks, native kernels,
jets, exponential boxes, or root searches. It independently audits the
new deductions and verifies saved arithmetic certificates.

## 1. A same-time event really gives a global positive floor

The 057 fixed basis is N=(J6-2I6)D with D=diag(-1,1,1,1,1,1). The stated
columns are body32 columns and detN=128. Since J6^2=6J6,

    (J6-2I6)(-I6/2+J6/8)=I6.

Thus N^-1=D(-I6/2+J6/8). Its diagonal absolute entry is3/8 and its five
off-diagonal absolute entries are1/8, giving row absolute sum1. This is
an infinity-norm bound; no Euclidean operator norm has been substituted.

Write the remaining columns as M. In the positive fixed-basis formula,
the 26 remaining parameters w retain their actual product Cauchy(t) law,
and the solved six parameters are N^-1(x-Mw). On ||x||infinity<=1 and
|w_j|<=t, each coordinate of Mw has magnitude at most26t, so every solved
parameter has magnitude at most1+26t. The event has probability exactly
(integral_-t^t t/[pi(t^2+w^2)]dw)^26=2^-26 for each positive t. Consequently

    p_t(x)>=2^-26/128 * [t/(pi(t^2+(1+26t)^2))]^6.

The basis formula is an everywhere-defined positive convolution density;
it is not a conditional law on a zero-probability fiber. This bound holds
on the entire cube, including S2. G may be infinite on S2, which is compatible
with the lower bound. No cutoff regularity or off-S2 assumption is needed.

For 1/32<=t<=1/16, an independent factored inequality is

    1765/16 - [1/t+52+677t]
      = (1/16-t)(677t-16)/t >=0.

Here 677t-16>=165/32>0. This verifies the majorant on the whole time
interval, rather than at two sampled times. Its reciprocal gives the
Cauchy-factor lower bound16/(1765pi). Including the probability, basis
Jacobian, all six factors and time-window length gives

    2^-26 * 128^-1 * (16/(1765pi))^6 * (1/32)
      =1/[2^14(1765pi)^6].

Using pi<22/7 yields the strict global floor

    G(x)>g_L=7^6/[2^14*38830^6].

Its t-dependent event is inside the correct t-density formula throughout.
No Cauchy parameter has been moved to another probability time. This proves
a complete but very small floor, not the candidate g0=1/65000.

## 2. Exact all-time integral and complete-face upper endpoint

The reviewed 036 pointwise-in-time estimate is

    ||p_t(1,.)||_(6/5,F)<=A*t^(1/6)/(1+1024t^2), t>0.

The previous report's conditional final cost does not make this inequality
conditional on a global g certificate. The five-dimensional projection peak,
true positive-face L1 factor1/32, and interpolation coefficient A were
already independently certified. Its original complete dyadic A and pi
records are the inputs here. They are neither recomputed nor rebuilt from
printed values.

Near0 the time envelope is O(t^(1/6)); near infinity it is O(t^(-11/6)),
so both ends are integrable. Set s=32t and then s=z^6. Exactly,

    J:=integral_0^infinity t^(1/6)/(1+1024t^2)dt
      =32^(-7/6) integral_0^infinity 6z^6/(1+z^12)dz.

Here is an independent precise version of the residue calculation. Write
zeta=exp(i*pi/12); the six poles in the upper half-plane are zeta^(2k+1),
k=0,...,5. Each residue is (1/2)zeta^(-5(2k+1)). Their sum S/2 has

    S*(zeta^5-zeta^-5)=2,
    S=-i/sin(5pi/12).

The first identity is checked algebraically modulo Phi24(w)=w^8-w^4+1.
Since the rational integrand is even, the positive half-line integral is
one half of the real-line contour integral: (pi*i/2)S=pi/[2sin(5pi/12)].
The large semicircle contribution is O(R^-5), so there is no omitted arc.
There are no real poles. Since sin(5pi/12)=(sqrt6+sqrt2)/4, positive roots
u^6=2 and v^6=3 give

    J=pi/[32*u^2*(v^3+1)].

This establishes both the time scaling and the factor1/2 from the even
contour. No numerical time quadrature or temporal grid is used.

Minkowski on [a,T], then increasing positive time truncations and monotone
convergence, proves ||G(1,.)||q<=A*J. It includes the full face and its
singular neighborhoods; possible infinite values on a null set are not
discarded by hand. In particular I is finite.

All quantities in J are positive. Inserting A_U and pi_U in the numerator
and the strict original sixth-root lower endpoints u_L,v_L in the denominator
therefore gives the valid strict bound B=I^(5/6)<B_U. A saved rational
I_U satisfying I_U^5>=B_U^6 then gives I<I_U. The exact root brackets, not
decimals or two separately widened intervals, are the proof objects.
This is a new complete-face consequence of 036, without changing its
historical local-only result.

## 3. The point upper barrier and real-box lower integral

The g_U input is precisely 032 results.json#/points/0/check/G_interval/1.
The frozen original point is in the cube and has a positive certified
distance to the true S2. The accepted052 theorem gives finite continuity
of G in its neighborhood. For every eta>0 an inward positive-volume part
of a sufficiently small neighborhood has G<=G(point)+eta. Hence the
essential infimum over the six-dimensional cube satisfies
g<=G(point)<=g_U. Merely evaluating an arbitrary boundary point would not
justify this without the neighborhood/continuity argument.

Both 061 lower boxes concern exactly the original body32 band
B_(1/250,1/2). The final 061 audit seal, completed through the twelve-field
reference work, binds the original jets and box files by hash. It preserves
the earlier exact native-to-jet and whole-box tangent audit and closes
the same two cases' gradient gap. It does not certify arbitrary points or
windows. The current audit checks those file bindings and the original
lower values; it does not treat the preparation-only partial reference as
complete and does not rerun any native computation.

The boxes are clipped, genuine five-dimensional face boxes. Their fourth
free-coordinate intervals are strictly separated. Their exact volumes are
products of the actual clipped widths, not (2/1000)^5 and not a sample
weight. Let their certified positive band floors be l_i and volumes V_i.
Since full G dominates this positive band everywhere,

    I>=sum_(i=1)^2 V_i*l_i^(6/5).

A lower fifth-root endpoint for l_i^6 is a valid lower bound for l_i^(6/5),
so summing volume times those lower endpoints gives I_L. Uncovered face
regions contribute nonnegative mass and are set to zero only in this lower
estimate. This cannot be confused with a full-face average or upper bound.
The positive band is smooth across its S3 anchor; no unproved full-G
derivative at that anchor is used.

## 4. Structural model floor without restricting a.e. data to a face

It is not legitimate just to restrict a six-dimensional a.e. inequality
to the five-dimensional boundary. Use the exact homogeneity
G(r x)=r^-5 G(x), which follows from the common-time stable scaling and
positive time substitution. Parameterize the positive cube sector in
which coordinate1 is maximal by (r,y)->r(1,y), 0<r<1, y in F. Its volume
Jacobian is r^5. Since G>=g almost everywhere in the cube, Fubini implies
that for almost every y, for almost every r in(0,1),

    G(1,y)=r^5 G(r(1,y))>=r^5 g.

Taking a sequence of such r tending to1 gives G(1,y)>=g for almost every
y. Therefore I>=g^q because |F|=1. The baseline g_L>0 and point barrier
g_U<infinity ensure this involves a finite positive g.

For the accepted same-generator positive-part comparison, optimizing the
constant truncation height c in (0,g) gives cost
5*c^(-1/5)*I/(g-c), minimized at c=g/6, namely
U_model=6^(6/5)*I/g^(6/5). Thus I>=g^q implies U_model>=6^q. This is a
floor for the fixed-potential truncation model, not for the maximal
operator C6. The existing full L1/all-scale operator-transfer interface
is not reconstructed by this scalar audit.

## 5. Monotone endpoint transfer and exact meanings

All endpoints belong to this one fixed model, with 0<g_L<=g<=g_U and
0<I_L<=I<=I_U. Because U increases with I and decreases with g,

    max(6^q,6^q I_L/g_U^q) <= U_model <= 6^q I_U/g_L^q.

The structural term is represented numerically by its saved lower root
endpoint, so the reported rational interval remains outward. Cost roots
have targets6^6*I_endpoint^5/g_endpoint^6. This verifies both the fifth
root exponent and the numerator/denominator choices.

For actual-model surplus S=64(g/6)^q-I, monotonicity gives

    64(g_L/6)^q-I_U <= S <= 64(g_U/6)^q-I_L.

The first power is rounded down and the second up. The supplied outward
U interval crosses64 and S interval crosses0. Neither establishes true
model exclusion or a competitive operator upper bound.

The conditional plan uses the unproved value g0=1/65000 with this fixed
I_U. Its U and S intervals enclose that plan expression only, not an
assertion that the actual global g equals or exceeds g0. Finally the plan
function 64(g/6)^q-I_U increases with g. Even its upper value at the rigorous
barrier g_U is negative, proving that improvement of g alone cannot rescue
this particular coarse I_U. It does not exclude the true G, because actual
I may be much smaller. This distinction matches the new guidance's use
of model-cost intervals rather than failures of a single upper envelope.

## 6. Finite independent checks and stopping scope

Four configurations are registered: the new exact integral identities;
the same-time event normalization and factored t majorant; direct power
checks for ten saved fifth-root certificates and two original sixth-root
brackets; then same-window/source/volume and endpoint-transfer checks.
Each root is validated directly by lower^n<=target<upper^n. No root is
searched again, no old six-configuration producer is executed, and no
036 peak or 061 kernel/jet/box exponential is recomputed.

All source dependencies and original windows are explicit and immutable.
No substantive mathematical gap has been found in the new 091 deductions.
The final acceptance is conditional on the fixed checker finishing; any
failure is saved and stopped without retry. Parent retains final proof
responsibility. There is no authorization for this subagent to push.
