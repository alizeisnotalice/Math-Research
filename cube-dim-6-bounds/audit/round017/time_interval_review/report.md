# F017 independent audit: continuous time interval obstruction

**PASS for body32 on J=[223/5000,251/5000], and for every probability
mixture of normalized short potentials supported in J.** At the frozen
B015 point x, K_T(x)/T<1/4100 throughout J. Consequently the standard
cube-domination model cost is at least 1025/16=64.0625 for every such
mixture. This is not a lower bound for the actual cube maximal constant.
No all-time, unrestricted-mixture, other-frame, or new C6 conclusion follows.

One ancillary normalization statement needs an explicit hypothesis:
mixing the raw K_T against a positive measure nu requires
0<integral T dnu(T)<infinity before dividing by kernel mass. Finiteness
of nu alone is insufficient on an unbounded time set. The approved
compact-J probability-mixture conclusion is unaffected.

## 1. Product measures, linear pushforward, and all-space density comparison

For c_t(z)=t/[pi(t^2+z^2)] and s>=t,

    c_s(z)/c_t(z)=(s/t)(t^2+z^2)/(s^2+z^2)<=s/t.

Taking the product over all 32 independent body-direction parameters gives
a genuine inequality of nonnegative measures on R^32,

    product_j c_s(z_j) dz <= (s/t)^32 product_j c_t(z_j) dz.

Push both measures through the same linear body-frame map A. For every
Borel set E, apply the measure inequality to A^-1(E); no conditional
distribution on a fiber is used. The pushforward measures have continuous
densities p_s and p_t because the frame spans R^6 and its Fourier symbol
is coercive. The measure inequality first gives the density inequality
almost everywhere. If it failed at a point, continuity would give an open
set of failure with positive measure, a contradiction. Therefore

    p_s(x)<=(s/t)^32 p_t(x), s>=t, every x.

Interchanging s and t gives p_s(x)>=(s/t)^32 p_t(x) for 0<s<=t. This
second inequality is a lower bound at earlier times. It is exactly what
is needed when subtracting the interval [T,t0] from K(t0).

The exponent is the number of original Cauchy factors, 32. Solving six
basis coordinates in a coarea formula does not delete the other 26
densities; nor does it create a time-dependent determinant factor. The
pushforward proof avoids both potential normalization errors.

## 2. The upper envelopes cover an entire continuum of times

Let K(t0) in [L,U], p_t0(x) in [l,u], with 0<l<=u, and t0=6/125.
Integrating the earlier-time lower bound gives, for 0<T<=t0,

    K(T)=K(t0)-integral_T^t0 p_s(x) ds
        <=U-l t0/33 [1-(T/t0)^33].

Integrating the later-time upper bound gives, for T>=t0,

    K(T)<=U+u t0/33 [(T/t0)^33-1].

These inequalities hold throughout their indicated half-lines, although
their envelopes are useful only on a limited interval. After division by
T, write the two expressions as

    E_left(T)=A_left/T+B_left T^32,
    A_left=U-l t0/33, B_left=l/(33 t0^32),
    E_right(T)=A_right/T+B_right T^32,
    A_right=U-u t0/33, B_right=u/(33 t0^32).

The saved rational data prove A_right>0. Since u>=l, A_left>=A_right>0,
and both B coefficients are positive. For either envelope,

    E''(T)=2A/T^3+32*31*B*T^30>0.

A convex function on a compact interval is bounded above by the maximum
of its endpoint values. Thus E_left on [a,t0] and E_right on [t0,b]
are controlled by the three rational numbers E_left(a), U/t0, E_right(b).
There is no inference of convexity or unimodality of K(T)/T itself.

The independently reconstructed values, displayed only for orientation,
are

    a=223/5000: E_left(a)  = 0.00024367340944126832...,
    t0=6/125:  U/t0       = 0.00023228107620985218...,
    b=251/5000: E_right(b)= 0.00024295653734983247....

Each is strictly below 1/4100 by a positive exact rational residual in
`results.json`. Hence the uniform strict inequality q(T)<1/4100 holds
at every T in the closed interval J, including its endpoints. The
recorded negative q'(t0) is consistent but is not used to extrapolate
monotonicity or a unique maximum.

## 3. Arbitrary supported probability mixtures and essential minima

For a probability measure mu supported in J, define

    Q_mu(y)=integral_J K_T(y)/T dmu(T).

Tonelli gives integral Q_mu=1. Order-one scaling gives

    r^-6 [K_T(y/r)/T]=K_(rT)(y)/(rT), r>0.

After convolution with nonnegative f, the dilation indexed by r is thus
a convex average of the same Markov semigroup's prefix-time averages.
It is bounded by their common ergodic maximal function, for every r.
This is the already audited positive Markov weak-(1,1) interface; no
supremum is interchanged with convolution to invent a new kernel.

The three-endpoint proof supplies a common rational M<1/4100 with
q_x(T)<=M for every T in J. Therefore Q_mu(x)<=M, even for nonatomic or
singular probability measures mu.

To pass from this point value to an essential cube minimum, reuse the
independently audited F015 exact geometry at

    x=(100000000000000,99999998999999,99999997999988,
       90268521812170,53714088833376,157107163435)/10^14:

    dist(x,S2)^2
      =6020805515960314567115330493/40000000000000000000000000000>0.

Choose a compact spatial neighborhood E disjoint from S2. The previously
proved positive-parameter greedy-rank estimate bounds p_t uniformly on E
for 0<t<=1. Fixed positive-time densities are jointly continuous; the
small-time uniform bound makes integral_0^T p_t continuous jointly in
(y,T) on E x J by dominated convergence. Since J is bounded away from
zero, K_T(y)/T is jointly continuous and bounded there. Dominated
convergence in mu therefore makes Q_mu continuous near x. Equivalently,
joint uniform continuity on the compact product supplies a neighborhood
estimate uniform over all these probability measures.

Although x lies on the face x1=1, every sufficiently small neighborhood
of x intersects the cube interior in a set of positive six-dimensional
measure. Thus essinf_cube Q_mu<=Q_mu(x)<=M<1/4100. The standard
domination cost 1/(64 essinf_cube Q_mu) is consequently at least
4100/64=1025/16; if the essential infimum is zero, that cost is infinite.
This does not assert a lower bound for C6 itself or rule out other proof
methods.

For a finite positive measure mixing already normalized K_T/T, division
by its finite nonzero total measure gives the probability formulation.
For raw K_T mixtures, the normalization mass is integral T dnu(T), not
nu((0,infinity)); the corresponding probability measure is
T dnu(T)/integral T dnu(T). A finite, positive first time moment is
therefore required. It holds automatically for finite nonzero measures
supported in this compact positive J. This is the only normalization
qualification to the ancillary wording in source theorem section 1.

## 4. Actual independent replay and source boundaries

The checker in this directory reads the frozen B015 native output
`02.stdout.json`, its checkpoint, C017's source manifest and delivery
receipt, and F015's prior independent result. It decodes the native
dyadic K and p endpoints exactly, checks the full 967680-term supported
status, and matches the exact point, endpoints, hashes, and previously
audited off-S2 geometry. It imports or executes no producer code.

It reconstructs the two envelope endpoints using the separated
A/T+B*T^32 formula, unlike the producer's integrated bracket form, then
matches the saved exact rational answers. It reproduces exactly the same
11 named scalar checks and saves every positive residual. There are no
extra spatial tests, new kernel times, G values, root computations, or
parameter searches. The low-precision display fields are not proof inputs.

The actual run exited 0: PID **35967**, Python **3.15.0**, script time
**0.0009922920144163072 seconds**, tool chunk **37567c**. It used Fraction
integer arithmetic only, with zero kernel, Monte Carlo, MATLAB, or Sturm
calls. Reproduce with

    /Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python verify_saved.py

The local verifier has the mathematical interpretation above; file hashes
bind identities but do not replace the previously completed F014/F015
kernel and spatial audits. All input/source identities are in the result
and manifest. The seven C017 source files and every earlier F directory
remain unchanged; no stage, commit, or push was performed.

## 5. No extrapolation to all times or other frames

As T approaches zero, the proved A_left/T envelope diverges. As T tends
to infinity, the proved B_right*T^32 envelope diverges. Their validity on
the half-lines does not extend the useful q<1/4100 certificate beyond J.
The previous fixed-time failures at .032/.048/.064 do not supply missing
time intervals, a unique peak, or a common all-time point obstruction.

The source's exact rational-density/root-isolation discussion is an
unexecuted conditional route. This audit does not expand its denominator,
run Sturm isolation, evaluate proposed critical times, or certify the
missing global hypotheses. An unrestricted-support probability-mixture
obstruction would first need a true uniform all-T point bound. Local
uniform domination could then justify continuity for that mixture, but
domination alone supplies no such all-T bound.

The present result concerns body32 only. It makes no statement about
augmented or pair-orbit frames. The current C6 upper bound remains 64.
F017 is complete and frozen; no round020 source was executed or evaluated.
