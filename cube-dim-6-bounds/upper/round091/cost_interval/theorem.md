# 091: certified model-cost interval and an explicit zero-kernel baseline

The model throughout this note is the original six-dimensional body32,
V={(1,+/-1,...,+/-1)}, independent unit Cauchy factors at the same time,
G=integral_0^infinity p_tdt. Put q=6/5, F=[0,1]^5,
g=essinf_([-1,1]^6)G and I=integral_F G(1,y)^qdy. The established positive-
part source and all-scale/general-L1 comparison give the fixed-potential,
optimal-constant-truncation model cost U=6^q I/g^q. It is a sufficient C6
upper certificate when evaluated with certified g_L and I_U, not a C6 lower
bound when evaluated from below. Existing C6<=64 remains the best checked
operator upper bound. Nothing here resolves the target g>=1/65000.

## 1. A coarse global g_L, without evaluating any old point

Use the already proved057 six-column body basis
N=(J6-2I6)diag(-1,1,...,1), detN=128. Its inverse is
diag(-1,1,...,1)(-I6/2+J6/8), and each row has absolute sum1. Write the
remaining26 columns as M. For every ||x||infinity<=1 and |w_j|<=t,

    |(N^-1(x-Mw))_i|<=1+26t.

Each remaining Cauchy(t) variable falls in [-t,t] with probability1/2.
The positive fixed-basis density formula therefore proves, pointwise on the
entire cube (including the exceptional set),

    p_t(x)>=2^-26/128 * [t/(pi(t^2+(1+26t)^2))]^6.

On [1/32,1/16], the convex function 1/t+52+677t is bounded above by the
larger endpoint value 1765/16 (the other is3365/32). Hence the bracket is
at least16/(1765pi). Integrating over length1/32 gives

    G(x)>=1/[2^14(1765pi)^6] > g_L:=7^6/[2^14*38830^6].        (1)

Only pi<22/7 was used in the last strict inequality. The event restriction
depends on t inside the correct same-time density formula; no factor has
been changed to a different probability time. This is deliberately a very
weak explicit global baseline, not a useful improvement of the model cost.

## 2. Complete I_U by integrating an already certified all-time envelope

036 theorem(6), independently reviewed there, holds for every t>0:

    ||p_t(1,.)||_(6/5,F)<=A t^(1/6)/(1+1024t^2),
    A=(32R1/pi^6)^(1/6)pi^(-5/6), R1=67320683/3251298631680.

Reuse its full saved Arb interval for A and pi, not displayed decimals or
a rerun of the projection peaks. Put u=2^(1/6), v=3^(1/6), both positive.
An exact positive-time integral is

    integral_0^infinity t^(1/6)/(1+1024t^2)dt
       =pi/[32u^2(v^3+1)].                                   (2)

For completeness, scaling s=32t followed by s=z^6 reduces the integral to
32^(-7/6) times integral_0^infinity 6z^6/(1+z^12)dz. Close the even rational
integral in the upper half plane; its six simple poles are
exp(i(2k+1)pi/12), k=0,...,5. Their geometric residue sum gives
pi/(2sin(5pi/12)), and sin(5pi/12)=(sqrt6+sqrt2)/4. Simplifying with u^6=2
gives (2). The contour arc vanishes by degree decay. This is an exact
identity, not a numerical quadrature or a new temporal grid.

Apply Minkowski first on finite time intervals, then monotone convergence.
It follows that B=I^(5/6) is at most A times(2). With the original complete
endpoints A_U, pi_U, u_L, v_L, set

    B_U=A_U*pi_U/[32u_L^2(v_L^3+1)].

Then B<B_U, and any rational I_U with I_U^5>=B_U^6 is a valid complete-face
upper endpoint. This closes a coarse numeric I slot with zero new kernels.
It is not a competitive integral bound, and is not an improvement of C6.
The old036 report's local status is preserved; this is a new all-time
consequence of its inequality, not a rewrite of its historical output.

## 3. Lower endpoints from complete original certificates

Use the FULL032 point interval endpoint, results.json#/points/0/check/
G_interval/1, rather than its shortened G_rational_display_bounds. The
exact off-S2 certificate and052 continuity give g<=g_U, including for an
essential infimum on the cube: an inward positive-volume neighborhood of
the boundary point has values arbitrarily close to its finite value.

The two real061 whole-box certificates and the completed065 independent
twelve-field check refer to the same band F_aT=integral_(1/250)^(1/2)p_tdt
and the same body32. Their original face_box intervals are already clipped
to [0,1]^5. They are disjoint in the fourth face coordinate. If their saved
strict lower bounds are l1,l2 and their exact volumes are V1,V2, positivity
G>=F_aT and monotonicity give

    I>=V1*l1^(6/5)+V2*l2^(6/5)>=I_L.                          (3)

Only new power/volume arithmetic is performed; no point or gradient is
re-evaluated. This very small covered volume makes(3) a weak model lower
endpoint. It must not be interpreted as a full-face mean or sampled integral.

There is also the structural floor U>=6^(6/5): g is at most almost every
face value. This follows from homogeneity and off-S2 continuity, or the
radial cube-coordinate formula identifying the essential infimum with its
face infimum. Thus I>=g^q|F|=g^q. This floor is about the constant-truncation
model, not the operator C6.

## 4. Cost and surplus tables

The complete, same-model endpoints give

    max(6^q,6^q I_L/g_U^q)<=U_model<=6^q I_U/g_L^q.             (4)

For the actual model surplus S=64(g/6)^q-I,

    64(g_L/6)^q-I_U<=S<=64(g_U/6)^q-I_L.                      (5)

The verifier encloses every fifth root by exact rational inequalities.
It separately records the conditional plan with g0=1/65000 (unproved),
and the best possible plan retaining this fixed I_U but allowing g up to
the rigorous g_U. If even the latter has negative surplus, improving g
alone cannot rescue THIS coarse I certificate. That is not an exclusion
of the true potential, whose true I could be much smaller. The interval
in(4), rather than a one-sided allocated budget, decides whether the
fixed potential has actually been excluded. No new C6 bound is claimed.

All inputs are complete saved certificates with hashes. Different alpha,
frame, or time-average models are separate registry rows. Their g/I cannot
be substituted here. The general-L1/all-positive-scale comparison is the
existing same-generator interface; this note changes no operator convention.
