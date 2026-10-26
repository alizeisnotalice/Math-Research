# 074: one fixed long-time middle-union budget

This is the original equal-weight body32 Cauchy model, q=6/5. Let rho=1/64,
epsilon=1/1000, and D be the union of the31 unsigned tubes with the B mean
in (1/4,3/4]. We safely include the endpoint1/4 in upper estimates. This
note estimates only H=integral_epsilon^infinity p_tdt on D. The051 short
certificate is read-only; no previous scalar or body kernel is recomputed.

## 1. Preserve the actual time terms in the complementary event

Use exactly the051 exhaustive parameter partition in each tube: both
specified parameters a_ones,a_w exceed1/2,1/16 in modulus, or its actual
disjoint complement. Read h_v0,h_v1,h_w0,h_w1,g_k and D_big,k from full
saved dyadic records. These already include every signed line/pair and
the true block perturbation, with k=1,...,5 and m=6-k.

The043 exact Cauchy average, for z!=0, is

 J_t(z,b)=(1+t/gamma)/[|z|^2+(lambda+L^2)t^2+2L^2t gamma],
 gamma>=sqrt(lambda)t/L, L=|b|.

Consequently, also when L=0 by direct evaluation,

 J_t(z,b)<=(1+L/sqrt(lambda))/(|z|^2+lambda t^2).

This strengthens the old time-free bound without taking an invalid
infimum over the shifted residual. The same extension of the retained
positive Cauchy integral and missing-parameter quotient loss from051
therefore gives the new all-time upper functions

 K=(1000/729) Crem F, F=1+sqrt(6/928),
 P_k(t)=D_big,k/[1024 t^2(1/4+t^2)(1/256+t^2)],
 B_k(t)=K/(g_k^2+928t^2) *
     [1/((h_v0^2+992t^2)(h_v1^2+960t^2))
      +1/((h_w0^2+992t^2)(h_w1^2+960t^2))],
 F_k(t)=P_k(t)+B_k(t).

The factor1024 in P_k follows from the frozen D_big=64 R2/pi^6 and
the full coefficient R2/(16pi^6). Both branches are actual positive
densities dominated by p_t; their sum is p_t. Thus p_t<=F_k(t) on D_A.
Every F_k is positive and decreasing. No branch reflection symmetry is
used, and no new third-distance claim along the parameter fiber is made.

For all t>0,

 F_k(t)<=C_k t^-6,
 C_k=D_big,k/1024+2K/(992*960*928).

## 2. Full-time true positive-face marginal, not its t=0 simplification

The accepted049 marginal at fixed coordinates1,j is

 p2_t(1,z)=2048t^2/[pi^2(((1+z)^2+1024t^2)((1-z)^2+1024t^2))].

The remaining four positive-coordinate slice has exactly1/16 of this
mass, for the full density only. Let L_m,U_m be the saved051 coordinate
limits, including the singleton m=1 case. Then

 L_m(t)=(128t^2/pi^2) integral_(L_m)^(U_m)
                dz/[((1+z)^2+1024t^2)((1-z)^2+1024t^2)]

is an upper bound for integral_DA p_tdy. For a=32t>0 a primitive of the
reciprocal product is

 Q_a(z)=log(((z+1)^2+a^2)/((z-1)^2+a^2))/(8(1+a^2))
       +[atan((z+1)/a)+atan((z-1)/a)]/[4a(1+a^2)].

Every D point has at least one free coordinate in [15/64,49/64]. A positive
union over only these five physical coordinates, independent of the31
tubes and with no disjointness assumption, gives

 L_D(t)<=min(5 L_star(t), t/[pi(1+1024t^2)]),

where L_star uses exactly those rational endpoints in the same primitive.
This uses symmetry only for full p_t. The count is5, not31 or an extra
reflection factor.

The tube volume bound, with true face measure and the maximal-coordinate
cone fraction, is

 V_k=(pi^2/4)rho^4 sqrt(m/k).

It uses B-mean interval length1/2, twice062's high interval length; it is
an upper volume, not an assertion that the clipped tubes are disjoint.
The multiplicities are n_k=binomial(5,k-1)=(1,5,10,10,5).

## 3. Three valid instantaneous norm envelopes

By positivity, the union inequality and p_t^q<=F_k^(q-1)p_t on each tube,

 ||p_t||_(Lq(D)) <= min{
   [sum_k n_k V_k F_k(t)^q]^(1/q),
   [sum_k n_k F_k(t)^(q-1)L_m(t)]^(1/q),
   F_max(t)^(1/6) L_D(t)^(5/6) }.

The first two may overcount overlaps positively. The last pays one
whole-union mass and a uniform pointwise bound. No cancellation, independence
or disjointness between overlapping tubes is assumed. Minkowski then bounds
the requested long norm by the positive-time integral of any such envelope.

## 4. Fixed rational interval integration and infinite tail

There is no spatial grid or semigroup/kernel evaluation. Freeze the time
interval endpoints a0=epsilon, a_(i+1)=min((16/15)a_i,1/4) until T=1/4.
For each slab [a,b], put r=b/a. Monotonicity of the denominator gives

 F_k(t)<=F_k(a), L_m(t)<=r^2 L_m(a),
 L_star(t)<=r^2 L_star(a),
 t/[pi(1+1024t^2)]<=r a/[pi(1+1024a^2)].

Use these bounds in the three expressions above, choose a rigorously
bounded minimum, and multiply by b-a. This is an interval upper Riemann
certificate based on analytic inequalities, not unverified numerical
quadrature or point sampling of the norm. The fixed finite partition is
not refined after seeing a result. The rest is paid analytically:

 integral_T^infinity ||p_t||q dt
 <=[sum_k n_k V_k C_k^q]^(1/q)/(5T^5).

All slabs and the final tail cover [1/1000,infinity) exactly. The total
is a rigorous upper budget for the requested long potential. If this fixed
upper expression exceeds the062 sufficient remaining D+F allowance, that
is only failure of this expression; it is neither a lower bound on the
actual long norm nor an exclusion of the model. No threshold adjustment,
extra parameter scan or repeated run follows that outcome.

After this local budget,051 short must still be added on the same D by
Minkowski before placing its q power in059. Far and global g remain unpaid;
no local expression here is a C6 or Cn bound.
