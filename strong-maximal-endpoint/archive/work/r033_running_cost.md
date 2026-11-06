# r033: exact running-cost Bellman reduction for the r032 payoff

Date: 2026-11-06.  This is a structural reduction only; it does not prove
the sharp Bellman lower bound needed for the endpoint coefficient.

Let `r=s-X` and

`psi_N(r)=1-exp(-r)-r/N` for `r>0`.

Define the nonnegative terminal penalty on the whole real line by

`a_N(r)=(1-exp(-r)-r/N)_+`.

For `r<=0` this is zero automatically.  Besides the root zero, `psi_N`
has a unique positive root `q_N`, determined by

`1-exp(-q_N)=q_N/N`.

It is positive on `(0,q_N)` and negative above `q_N`; for large `N`,
`q_N in (N-1,N)` and `N-q_N=N exp(-q_N)` is exponentially small.

## Exact expectation identity and exponentially small correction

For the payoff from r032,

`H_{N,s}(X)=exp(-r)` if `r>0`, and `1-r/N` if `r<=0`,

so, pathwise,

`H_{N,s}(X)=1-r/N-psi_N(r)1_{r>0}`.

Let a policy choose `c_i` predictably, let `m_i=m(c_i)`, and put

`D_N=sum_i(mu-m_i)`.  Since `E X=N mu-E D_N`,

`E_pi H_{N,s}(X)
 =1+(N mu-s)/N-E_pi D_N/N-E_pi[a_N(s-X)] + e_pi,N(s)`,             (1)

where

`e_pi,N(s)=E_pi[(-psi_N(s-X))_+ 1_{s-X>0}]>=0`.

For `s=N mu-d_N+O(log N)`, this correction is uniformly exponentially
small.  Indeed `s-X>q_N` implies

`X<s-q_N <=(mu-1)N-d_N+O(log N)+1 <=alpha N`

for a fixed `alpha in (mu-1,1)` and all large `N`.  Couple every increment
as `Z_{c_i}=Y_i+c_iB_i` with `Y_i~Exp(1)`.  Then pathwise
`X>=G_N=sum Y_i`; hence the Gamma Chernoff bound gives

`P_pi(s-X>q_N)<=P(G_N<=alpha N)<=exp(-cN)`.

On this event `0<=(-psi_N(s-X))_+<=s/N+1=O(1)`, proving
`0<=e_pi,N(s)<=C exp(-cN)` uniformly in the policy.

## Equivalent min-Bellman problem

Define, for `t` remaining steps,

`V_0(r)=a_N(r)`,

`V_t(r)=inf_{c in [0,infty]}
 { (mu-m(c))/N + E[V_{t-1}(r-Z_c)] }`.                            (2)

Predictable dynamic programming gives

`V_N(s)=inf_pi E_pi[ D_N/N+a_N(s-X)]`.

Consequently (1) gives the uniform equivalence

`sup_pi E_pi H_{N,s}(X)
 =1+(N mu-s)/N-V_N(s)+O(exp(-cN)).`                               (3)

The use of a *minimum* Bellman problem is essential: high payoff corresponds
to small accumulated mean-loss plus terminal penalty.

## Threatening-policy localization

At the prospective sharp threshold

`s=s_Q+(g+epsilon)log N`, `d=N mu-s`,

the r032 desired all-control payoff estimate is equivalent, by (3), to

`V_N(s)>=d/N+c_epsilon(log N)^2/N+o((log N)^2/N)`.                (4)

No such bound is proved here.  Nevertheless (3) localizes the only policies
that can threaten an endpoint upper bound: a policy with

`E[D_N/N+a_N(s-X)] > d/N+O(log N/N)`

already has `E H_{N,s}(X)<1-O(log N/N)` and cannot support an input with
`E h(X)>1` after the r032 profile reduction.  Analysis may therefore focus
on near-minimizers of (2), where running drift loss and the terminal
shortfall penalty must balance.  This statement neither proves (4) nor a
sharp upper bound for `R_N`.

## Appendix: exact front identity (conditional numerical diagnostic)

Let `Gv(r)=int_0^infty e^-y v(r-y)dy` and write `A=Gv`.  For the fixed-
`c` averaging operator,

`K_c v(r)=[A(r)+e^-c A(r-c)]/(1+e^-c)`.

At an **interior** optimizing action `c_*`, differentiation and
`A'=v-A` give

`partial_c K_c v(r)=e^-c/(1+e^-c)[K_c v(r)-v(r-c)]`.

Thus the first-order condition is exactly

`T v(r)=v(r-c_*)`.                                               (F)

If `v` is strictly decreasing, `v(q_prev)=1`, `T v(q_new)=1`, and this
interior optimizer is unique, then

`c_*=q_new-q_prev`.

Substitution into the average identity yields

`e^qnew(A(qnew)-1)=-e^qprev(A(qprev)-1)`,

and, since `(e^u(A(u)-1))'=e^u(v(u)-1)`, equivalently

`int_qprev^qnew e^u(1-v(u))du
 =2e^qprev(A(qprev)-1)`.

For the initial terminal payoff from r032, `q_prev=0` and direct
integration gives

`A_0(r)=e^-r(r+1+1/N)` for `r>=0`.

Hence the first positive front is the unique positive solution of

`e^q1=1+q1+2/N`.

Under the stated Bellman regularity, this is in fact a global identity, not
only a diagnostic.  Define the first-order numerator

`D(c)=A(r)+e^-c A(r-c)-(1+e^-c)v(r-c)`.

If `v` is locally absolutely continuous and strictly decreasing, then,
using `A'=v-A`,

`D'(c)=(1+e^-c)v'(r-c)<=0` a.e.

and the integral of this derivative is strictly negative on every nontrivial
interval.  Thus `D` is strictly decreasing.  Also
`D(0)=2(A(r)-v(r))>0`.

For `r>=0`, take `c=r+mu`, so `r-c=-mu`.  The negative-ray continuation
`v(x)=1-x/N+(k-1)mu/N` gives
`A(-mu)=v(-mu)+1/N` and
`A(0)=v(-mu)+(1-mu)/N`.  Since `A` is decreasing,

`D(r+mu)<= [1-mu+e^{-(r+mu)}]/N<=0`,

where the final inequality is `mu=1+e^-mu` and `r>=0`.  Hence there is a
unique global maximizer, it is interior except for an equality endpoint, and
`0<c_*<=r+mu`.  The sign identity
`partial_cK_cv=e^-c D(c)/(1+e^-c)^2` proves global increase then decrease.

For the actual Bellman iterates, the terminal payoff has a positive local
decrease modulus on every bounded-above `r` interval:
`phi(r1)-phi(r2)>= (r2-r1)min(1/N,e^-max(r2,0))` for `r1<r2` and all
nonnegative increments.  The Bellman supremum preserves this local strict
decrease and the terminal Lipschitz bound preserves local absolute
continuity.  On the negative ray it is exactly the displayed affine
continuation.  Therefore the hypotheses used above hold for these iterates.
