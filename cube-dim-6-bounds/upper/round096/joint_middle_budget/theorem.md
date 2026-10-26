# 096: sharp time–spatial-rank coupling, then a finite information test

The model is the unchanged original body32 Cauchy semigroup, q=6/5,
Lebesgue positive face, and the31 middle tubes D_A at rho=1/64 and
mean_B in(1/4,3/4]. All time is t>0. Actual overlap is allowed. This
round introduces a different functional from074's integral of instantaneous
norm bounds. No old time table is refined or replayed.

## 1. A self-contained sharp rearrangement bound

Let X have measure at most V, padded by zero to a non-atomic space of
measure V. Suppose a nonnegative family f_t satisfies, for every t>0,

  0<=f_t<=F(t),    integral_X f_t <= M(t),

where F(t)>0. Put theta(t)=min(V,M(t)/F(t)) and

  h(s)=integral_0^infinity F(t) 1_{s<theta(t)} dt,  0<s<V.

Then for q>1,

  integral_X (integral f_t dt)^q <= integral_0^V h(s)^q ds.       (R)

Proof: for any measurable A of measure a, Tonelli and the two bounds give

  integral_A integral f_tdt <= integral min(aF(t),M(t))dt
                            = integral_0^a h(s)ds.

The right equality uses a<=V. The function h decreases. For lambda>=0,
the identity integral(H-lambda)_+=sup_A[integral_A H-lambda|A|]
and the preceding bound imply integral(H-lambda)_+ <= integral(h-lambda)_+.
Integrate this inequality with weight q(q-1)lambda^(q-2), using the
elementary identity z^q=q(q-1) integral_0^z (z-lambda)lambda^(q-2)dlambda.
Positive truncations and monotone convergence justify possibly infinite
quantities. This proves(R) without differentiating events or assuming
any spatial/time independence.

It is sharp for this abstract information class: on(0,V), take
f_t^*(s)=F(t)1_{s<theta(t)}. This has mass min(VF,M), obeys the cap,
and its potential is h. Sharpness is on the abstract space, not a claim
that a Cauchy density can attain it.

Compared with interpolating at each time and applying Minkowski, (R)
retains the changing occupied fraction theta(t)/V inside the common
spatial rank integral. Minkowski applied to f_t^* proves that (R)'s
q-root is no larger than its own one-branch instantaneous norm integral.
It need not dominate074's time-varying minimum across different branches.

## 2. Two applications retaining both tube and whole-union information

Use exactly the positive F_k(t), the true two-coordinate L_m(t), and
V_k from074. Their underlying geometric/event proof is051/043 and049.
No density or scalar source is reevaluated.

For one tube define h_k by (F_k,L_m,V_k) and
J_k=integral_0^V_k h_k(s)^q ds. Permutation multiplicities are
n_k=binom(5,k-1)=(1,5,10,10,5). Positive union gives

  integral_D G^q <= J_tube := sum_k n_k J_k.                    (T)

For the whole D also use

 F_*(t)=max_k F_k(t), V_*=sum_k n_k V_k,
 M_*(t)=min(5L_star(t),t/[pi(1+1024t²)]),

with the exact positive-face L_star projection interval[15/64,49/64].
The actual union has volume at most V_*, so (R) gives a second bound
J_whole. Thus the new candidate is

  integral_D G^q <= J_new := min(J_tube,J_whole).               (N)

The common-time dependence of each F/M pair is kept inside h before
taking q. This is a new spatial-rank/time coupling, not the excluded
074 expression E(t). The two relaxations are still solved separately;
this is not the optimum under all simultaneous tube-overlap constraints.

The formulas cover all time and are finite. At t->0, each F is bounded
above by C(t^-2+1) and below by d t^-2 for positive constants, while
each M is bounded above by C' t². Hence theta<=C''t^4 and h(s) is
O(1+s^-1/4) as s->0. Its q power is integrable since q/4=3/10<1.
At infinity F=O(t^-6), which is integrable. The same argument applies
to max F_k and the whole mass. Finite existence is not a competitive
numerical certificate.

## 3. A fixed-window feasible witness for each sharp relaxation

Reuse only the frozen092 rational lower coefficients f_k,v_k,c_J,c_W.
Their already proved scope is0<t<=b=1/64, and they satisfy

 F_k(t)>=f_k>0, V_k>=v_k>0, L_m(t)>=c_J t²,
 F_*(t)>=max f_k, V_*>=sum n_k v_k, M_*(t)>=c_W t².

092's lower coefficients were proved on0<t<=b even though its final
lower integral began at1/1000. No new endpoint or old function is
computed here.

For any one of the six triples(f,v,c) just described, define on(0,V)

  fhat_t(s)= f 1_{0<t<=b} 1_{0<s<min(v,c t²/f)}.

Its height is <=F, its mass is <=c t²<=M on the window, and it is zero
at all other times. It is an admissible family in the precise abstract
relaxation. Its potential is

  Hhat(s)=f(b-sqrt(fs/c))_+ 1_{0<s<v}.

Consequently the sharp relaxation's objective is at least

  W(f,v,c)=2c f^(q-1) integral_0^m y(b-y)^q dy,
  m=min(b,sqrt(vf/c)).                                       (W)

The exact primitive yields

 W=2c f^(q-1) [ b^(q+2)/((q+1)(q+2))
               - b(b-m)^(q+1)/(q+1)+(b-m)^(q+2)/(q+2) ].

For lower arithmetic replace m by a certified rational lower endpoint;
the integral increases with m. Bound the positive powers outward and
respect the minus sign. No full F_k(b), L_m(b), old grid, or native
kernel is evaluated.

Let W_k^- and W_*^- be these certified witness lower values. Then

 J_new >= min(sum n_k W_k^-, W_*^-)=:W_new^-.                 (L)

Both sides of the outer min have been bounded from below. A lower
bound for only one candidate is insufficient.

## 4. What the certificate can reject

Use the exact062 saved D+F allowance and I_paid_CEH. With far optimistically
zero and the unproved conditional g0=1/65000, this proposal would cost

 I_plan=I_paid_CEH+J_new,
 U_plan=6^(6/5) I_plan/g0^(6/5),
 S_plan=64(g0/6)^(6/5)-I_plan.

There is no extra051 short term here: (N) covers the full G at all t>0.
If W_new^- exceeds the remaining I allowance (or its B equivalent),
and the direct certified S upper is negative, this fixed pair of sharp
relaxations cannot meet that allocation even with exact integration.
The witness is not a semigroup, not a body32 density, and not a C6
lower-bound example. It excludes neither the actual model nor use of
stronger joint constraints, actual tube volumes/overlap, extra marginals,
or a different g/I pairing.

If the witness does not separate the threshold, report UNRESOLVED and
stop. No extra window, parameter scan, or adaptive refinement follows.
