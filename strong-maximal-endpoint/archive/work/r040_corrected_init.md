# r040: corrected right-log uplift — fixed-time initialization audit

Date: 2026-11-06.  Define for positive saddle parameter

`F(eta)=I_mu(y(eta))-I_*(y(eta))-D(eta)`,

and extend it by zero for `r<=t mu`.  Adding `t F(eta(r/t))` to the
logarithmic uplift changes the positive-side rate from
`I_mu-alpha D` to `I_*-(alpha-1)D`, which is strictly below the coarse
positive rate for every `eta in (0,1)` when `alpha>1`.

## Near the centre

The r039 coefficients give

`F(eta)=K4 eta^4+O(eta^5)`.

Since `eta` is a smooth multiple of `y-mu` near the mean, the zero extension
is `C^3`.  It does not change the small-r or left-affine initialization
proof.  On a fixed compact positive `r` interval its value is bounded; a
larger, still fixed, density coefficient `s0` can compensate its negative
part because `f_k` has a uniform positive minimum there.

## Far-right asymptotic needed for the fixed-k proof

For the fixed reference action,

`M_mu(eta)=(1-eta)^(-1)(1-p_mu+p_mu e^(eta mu))`

has coefficient `2/mu` as `eta` tends to 1.  The maximizing factor in
`M_*` is

`(1+exp(-(1-eta)c))/(1+exp(-c))`.

It is at most 2 and, on choosing `c=(1-eta)^(-1/2)`, tends to 2. Hence

`M_*(eta)=2(1-eta)^(-1)(1+o(1))`,
`M_mu(eta)=(2/mu)(1-eta)^(-1)(1+o(1))`.

The Legendre saddle therefore gives

`I_mu(y)-I_*(y)=log(mu)+o(1)`,
`D(eta)/y -> p_mu/8`,
`F(eta)/y -> -p_mu/8`.                                             (1)

This can be made strict by sandwiching the Legendre transform: the
upper/lower bounds on `M_*` force its maximizer to `eta->1`, then give the
coefficient 2.  No unproved assertion about an exact finite optimizer is
needed.

For fixed `k`, r039's tail identity gives

`J_k/(r v_k) -> q_k=(p_mu/8)(k+H_k)/(k+1)`.

Thus the corrected logarithmic inner tail has enhancement

`alpha q_k-p_mu/8>0` for `alpha>=1,k>=4`; its minimum is already positive
at `k=4`.  Local Wolfram evaluation gave `p_mu=0.2178117057` and
`q_4-p_mu/8=0.0058990670`; it also verified the optimizing Bernoulli factor
tends to 2.  By continuity the same is uniform for all sufficiently large
matched `a_N`.

Choose a fixed `0<kappa<min(1,alpha q_k-p_mu/8)`.  The corrected candidate
is bounded below by `c r^k exp(-(1-kappa)r)`, whereas the old uniform
screen `phi_N(u)<=exp(-theta u)` with
`theta=1-kappa/2` bounds every feedback payoff by
`M_*(theta)^k exp(-theta r)`.  Hence the right tail remains dominated.

## Scope

This validates the proposed fixed-time initialization repair, including
uniformity in large matched `N`.  It does not prove a one-step global
supersolution in the intermediate time range; that would require new
control of the correction under every `K_c`.
