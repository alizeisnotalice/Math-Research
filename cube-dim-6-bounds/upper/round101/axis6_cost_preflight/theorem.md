# 101: an exact failure certificate for the six-axis Green-positive-part family

## 1. Fixed real model and the exact cube floor

Let c_t(s)=t/[pi(t^2+s^2)] and p_t(x)=prod_(i=1)^6 c_t(x_i), all factors
at the same positive time t. This is the actual product-Cauchy probability
semigroup with L=-sum_i|D_i| and symbol psi=sum_i|xi_i|. Its kernel is even,
coordinatewise decreasing in |x_i|, and p_t<=pi^-6 t^-6. Its exact scaling is
p_t(rx)=r^-6 p_(t/r)(x), hence G(rx)=r^-5G(x) for G=int_0^infinity p_tdt.

The product makes the cube floor transparent: for every |x_i|<=1 and t>0,
p_t(x)>=p_t(1,...,1). Integrate positively, allowing +infinity on singular
strata. Thus the minimum is at all corners, and

    g=inf_cube G=essinf_cube G
     =pi^-6 integral_0^infinity t^6/(1+t^2)^6dt
     =pi^-6 * (1/2)B(7/2,5/2)=3/(512pi^5).

The Beta value already occurs in the historical axis pilot; it is not a new
kernel value here. The new logical step is its exact global-minimum status.
For completeness Gamma(7/2)=15sqrt(pi)/8, Gamma(5/2)=3sqrt(pi)/4 and
Gamma(6)=120 give the displayed coefficient. Near a corner all six coordinates
stay away from0; small time is O(t^6), large time O(t^-6), so G is continuous.
Positive-volume cube neighborhoods approaching the corner show that the
essential infimum is the same value, rather than merely a boundary minimum.

At a point with k nonzero coordinates, small-time p_t is a positive multiple
of t^(2k-6). Consequently the true infinite Green set is support size<=2;
it is not removed from the distributional or integral proofs below.

## 2. Full-time whole-face norm and a finite true cost upper

Put q=6/5, I=int_[0,1]^5 G(1,y)^qdy and B=I^(5/6). Minkowski on positive
time truncations followed by monotone convergence gives

    B<=int_0^infinity c_t(1)[int_0^1 c_t(s)^qds]^(5/q)dt.

For0<t<=1, the actual one-coordinate integral is bounded by

    int_0^1 c_t^q <=||c_t||infinity^(q-1)*int_0^infinity c_t
                   =(1/2)(pi*t)^(-1/5).

Here the positive half-line mass is1/2, with no extra orthant factor. Also
c_t(1)<=t/pi. Therefore the short-time norm integrand is at most
2^(-25/6)pi^(-11/6)t^(1/6), integrable at0, including every face singular
layer and intersection. For t>=1, p_t(1,y)<=pi^-6t^-6 pointwise on a face of
volume1. It follows that

    B<= (6/7)2^(-25/6)pi^(-11/6)+1/(5pi^6).

Use the certified elementary bounds pi>157/50>3 and pi<22/7. Since
(157/50)^11>8^6, the first fractional pi power is <1/8; also
2^(-25/6)<1/16. Thus, strictly,

    B < B_R:=3/448+1/3645 <1/128.

This is an actual whole-face upper bound, not a norm estimate at sampled
points. Together with the exact g it implies

    U=(6B/g)^(6/5) < [8(22/7)^5]^(6/5) <12000.

All power comparisons used below reduce to integer powers of positive
rational numbers. This very loose upper endpoint exists only to close a
complete cost interval; no additional integration is warranted once the
strict true-cost lower endpoint below exceeds64.

## 3. A fixed-window true mean forces the actual cost above64

The face has volume1, so Jensen gives B>=A:=int_[0,1]^5G(1,y)dy. By positive
Tonelli and the exact Cauchy marginal,

    A=int_0^infinity [t/(pi(1+t^2))]
                      *[atan(1/t)/pi]^5dt.

On the fixed window t in[1/4,1/2], 1/t>=2>sqrt(3), so atan(1/t)>pi/3.
The entire integrand is nonnegative outside the window. Hence

    A > (1/243)int_(1/4)^(1/2) t/[pi(1+t^2)]dt
      =log(20/17)/(486pi).

For x=3/17, log(1+x)>x-x^2/2=93/578. This follows by integrating the
positive derivative residual x^2/(1+x); it is not a fitted logarithm value.
Canceling the exact pi factors against g gives

    6B/g >=6A/g
      > [512/243]*(157/50)^4*(93/578)=:R_L.

The registered exact rational comparison R_L^6>66^5 proves

    U=(6B/g)^(6/5)>66>64.

In particular this lower endpoint applies to the true integral I, not to an
upper-budget formula. The inherited old Sobol value is not used in any
inequality. A coarse explicit true-I interval may also be given:

    1/100000 < I <1/300,

by comparing [93/(578*486*(22/7))]^6 with (1/100000)^5 and (1/128)^6
with (1/300)^5. The decisive result is the much more informative correlated
cost interval66<U<12000, since g is exactly known.

## 4. Full-space source and centered-cube/general-L1 transfer

The six-dimensional finite-spanning-Cauchy proof in001 is rechecked with the
axes constants, rather than importing the body32 value1/(160pi^6):
psi_axis is between |xi| and sqrt(6)|xi|, Dom_2(L)=H1, and C_axis=pi^-6 in
p_t<=C_axis*t^-6. Thus the same transient split gives G locally L1 and
R_c=(G-c)_+ globally L1 for every c>0. Alternatively the finite I proved in
section2 and the positive cube-radial computation below give R_c in L1
directly, including unbounded thin channels along singular axes/pairs.

For test functions, the genuine coordinate Cauchy second-difference formula
gives Lphi in L1 intersect L-infinity intersect C0,
||Lchi_R||infinity=O(R^-1), ||Lchi_epsilon||1=O(epsilon^5), and
(P_tphi-phi)/t ->Lphi in L-infinity. Pairing G against |Lphi| is legitimate
because int_0^infinity min(||Lphi||infinity,C_axis*t^-6||Lphi||1)dt<infinity.
It follows on the WHOLE space that LG=-delta0.

Positive semigroup Tonelli gives G=P_tG+U_t a.e., U_t=int_0^tp_sds. The
convex positive part and its1-Lipschitz property give P_tR_c>=R_c-U_t.
Test, divide by t and pass to0 to obtain the positive distribution

    mu_c=LR_c+delta0>=0.

It is a positive Radon measure. Large cutoffs and R_c in L1 give total mass1.
Writing h=min(G,c), mu_c=-Lh; the epsilon^5 small-cutoff estimate and h<=c
give mu_c({0})=0. This includes every singular plane and its overlaps; no
pointwise claim on the singular-set complement has been extended by fiat.

For clarity, the obstacle input from001 applies here with C_axis: for bounded
compactly supported f>=0, transient common-time semigroup approximation gives

    f=sigma-Lu, u>=0 in L1 intersect L-infinity intersect H1,
    0<=sigma<=kappa, sigma=kappa on D={u>0}, int sigma=int f.

Its proof only needs the actual t^-6 bound, spectral comparability, positivity
and conservation, all just checked for axes6. The discrete-potential bound is
F+(6/5)(C_axis*M)^(1/6)F^(5/6), not the body32 numeric constant. The same
spectral high-frequency tail and large-cutoff mass recovery apply. No theorem
about six independent rectangular radii is involved.

For0<c<g, K_c=R_c/[64(g-c)] is even, nonnegative, L1 and >=1_cube/64.
Order-one scaling and the full-space finite-source identity give

    K_(c,r)*Lu=(mu_(c,r)*u-u)/[64(g-c)r],
    K_(c,r)(x)=r^-6 K_c(x/r).

On the SAME D complement this is nonnegative. Therefore K_(c,r)*f is bounded
there by kappa*int K_c for every fixed r. Merge the positive-rational-r null
sets, use L1 dilation continuity and bounded f to cover all real r>0, then
monotone-truncate |f| to extend to general L1. The cube average is dominated
at exactly this one scalar r. The operator remains centered, axis-parallel,
common-radius and all-scale; a product semigroup has not been replaced by a
composition of one-dimensional maximal operators or a strong maximal operator.

Cube-radial coordinates have radial Jacobian r^5. Coordinate signs and the
six maximal-coordinate faces give J=384I (ties null). Direct positive radial
integration yields

    int R_c=384I*(5/6)c^(-1/5)=320I*c^(-1/5),
    int K_c=5I*c^(-1/5)/(g-c).

The logarithmic derivative is -1/(5c)+1/(g-c); the unique minimum is c=g/6,
with both endpoint costs diverging. Consequently the ideal family cost is
exactly U=6^(6/5)I/g^(6/5). Its strict lower bound above66 excludes every
spatial truncation c of THIS equal-axis Green-positive-part family from
improving64, even with exact I and exact g. It is not a lower bound on C6,
not an exclusion of arbitrary other kernels/corrections from the axes
semigroup, not an anisotropic weight-family exclusion, and not body32 work.

## 5. Fixed run and stop rule

One registered pure-exact program checks source handoff/cost normalizations,
the fixed-window strict lower, and the full-time upper. No old corner kernel,
old Sobol/time quadrature, body native, J100, MATLAB or recursive agent is run.
The already-proved E098 rational pi certificate is read-only; no old point
evaluation is hidden as a new configuration. No parameters, window, grid or
precision are expanded after the strict failure. Once the registered interval
is verified, seal this limited noncompetitiveness result and stop.
