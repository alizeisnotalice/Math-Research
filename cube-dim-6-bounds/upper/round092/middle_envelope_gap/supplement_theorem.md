# 092 authorized replacement: one analytic lower window

The original requested endpoint reuse cannot be carried out: frozen074 saves
only the selected upper rectangle interval on each of86 slabs. It saves no
F_k, L_m, L_star or three-candidate endpoint arrays. The selected candidate's
lower endpoint need not bound the minimum of the unsaved other candidates;
074 also replaces F_max by an exact rational upper endpoint. No endpoint is
reconstructed from a decimal or by replaying074.

Root explicitly authorized this replacement: fixed b=1/64, at most3 new
configurations, CPU30/program wall180 including cleanup, single attempt.
All previous074 sources stay read-only. This proves a lower bound on the
same positive mathematical envelope E, never on the true density norm.

Write the three074 branches as

 E_V(t)=(sum n_k V_k F_k(t)^(6/5))^(5/6),
 E_J(t)=(sum n_k F_k(t)^(1/5)L_m(t))^(5/6),
 E_W(t)=(max F_k(t))^(1/6)
         min(5 L_star(t),t/[pi(1+1024t²)])^(5/6),
 E(t)=min(E_V,E_J,E_W).

## Shared positive-coordinate integral lower bound

Saved exact endpoints certify [1/4,3/4] is a subset of every [L_m,U_m],
including equality for m=1. It is also a subset of [15/64,49/64].
For z in that fixed interval and 0<t<=b=1/64, let u=z² and c=1024t²<=1/4.
The denominator product is

 D(z,t)=(1+u+c)²-4u.

It increases with c; its u derivative is 2(u+c-1)<=-3/8 on this rectangle.
Consequently its maximum is attained at z=1/4,t=b and equals377/256.
By integrating the nonnegative integrand over the common half-length
interval, all five L_m and L_star satisfy

 L(t)>=(16384/(377 pi²)) t² = c_J t².

Also t/[pi(1+1024t²)] >= t²/[b pi(1+1024b²)]
= (256/(5pi))t². Thus the whole-union mass expression is at least
c_W t², where c_W=min(5c_J,256/(5pi)). Both branches of this inner min
are bounded from below before taking their minimum.

## Uniform pointwise-formula lower bounds; no old function replay

The F_k formula is P_k+B_k with both positive. Drop P_k entirely. For each
k use the saved upper endpoints of h_v0,h_v1,h_w0,h_w1,g, and the saved
lower endpoints of Crem and the average factor to form K_-.
All denominator factors increase with t, hence for 0<t<=b

 F_k(t) >= f_k := K_-/(g_+²+928 b²)
   * {1/[(h_v0,+²+992b²)(h_v1,+²+960b²)]
      +1/[(h_w0,+²+992b²)(h_w1,+²+960b²)]}.

This is a new one-window lower coefficient, not evaluation or refinement
of074's86 old endpoint functions. Downward rational truncation of f_k
is allowed and recorded. The saved lower endpoints of V_k are used.
Every power/root below is enclosed using pure integer inequalities.

Define positive rational lower coefficients A,J,W for

 (sum n_k V_k f_k^(6/5))^(5/6),
 (c_J sum n_k f_k^(1/5))^(5/6),
 (c_W (max f_k)^(1/5))^(5/6),

respectively. Then each of the three branches has its own lower bound,
so with C=min(J,W),

 E(t) >= min(A,C t^(5/3)) on [epsilon,b], epsilon=1/1000.

Outside this window use0. This gives a legal lower integral over the entire
[epsilon,infinity), while the unchanged074 total remains an upper bound.
No assumption about which074 branch won a slab is used.

## Exact integration and correct cost notation

Let x=(A/C)^(3/5). If epsilon<x<b, the exact lower-envelope integral is

 J0 = A*b - (5/8)*A*x - (3/8)*C*epsilon^(8/3).

The other two cases are A(b-epsilon) for x<=epsilon and
(3/8)C(b^(8/3)-epsilon^(8/3)) for x>=b. Here epsilon^(8/3)=1/10^8 and
b^(8/3)=1/65536 exactly. Outward root endpoints and the signs in this
formula give a rigorous [J0_lower,J0_upper]. Ambiguous switch placement
or failed arithmetic terminates the one attempt.

Set L=int_epsilon^infinity E(t)dt, R=old074 certified rectangle-plus-tail
upper, and c=old051 short norm certificate. Thus J0_lower<=L<=R.
If J0_lower exceeds the saved062 allowance, no exact time integration of
this same E can satisfy that sufficient allocation. This alone is not a
model rejection. The stronger direct conditional budget comparison is

 I_plan(L)=I_paid_CEH+(c+L)^(6/5),
 U_plan=6^(6/5)I_plan/g0^(6/5),
 S_plan=64(g0/6)^(6/5)-I_plan,     g0=1/65000 (unproved globally).

The unknown far contribution is set to0 only for this optimistic plan
diagnosis, not declared paid. A negative upper endpoint of S_plan rejects
this fixed conservative accounting with E even in that optimistic case.
It does not lower-bound actual I, actual model cost, R_c, or C6.

The maximum recoverable long-budget reduction from perfect integration is
at most R-J0_lower. The possible improvement in S from replacing the old
rectangle plan is at most I_old_plan,upper-I_plan,lower. These are upper
bounds on recoverable error, not an exact measurement of quadrature loss.

Notation correction: the first saved-cost configuration was run before the
full guidance was supplied. Its fields named S actually contain U. Preserve
the receipt and numerical data; the supplemental results use the guidance's
correct S margin and explicitly identify that historical field-label error.
