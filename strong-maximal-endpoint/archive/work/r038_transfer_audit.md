# r038: late initialization area cannot by itself smooth to reference scale

Date: 2026-11-06.  Let `T=N^gamma`, `n=N-T`, and let the endpoint be
`q=N mu-d_N`, with `d_N asymp sqrt(N log N)`.  The r037 late splice gives a
global one-step barrier after `T`, but an initialization error must still be
propagated by the controlled MAX semigroup.

## Explicit obstruction to an L1-only transfer

Even require the error to be bounded and supported inside the natural tube
at time `T`.  Put `y0=T m(a_N)` and, for `M=c log T`, take for example

`E_T(y)=e0 1_[y0,y0+M/e0](y)`,

with fixed small `e0>0`.  It has `0<=E_T<=e0`, mass `M`, and its support is
inside the `T`-tube for large `N`.

Choose a *constant* future action `c_N` such that

`m(c_N)= (q-y0)/n = mu-(d_N-T delta)/n`,

where `delta=mu-m(a_N)=O(log N/N)`.  This action exists near `mu`: since
the strict maximum of `m` is quadratic,
`|c_N-mu|=O((log N/N)^(1/4))`.  The Gamma--Bernoulli local limit theorem is
uniform on that compact parameter neighbourhood, so its `n`-fold density
at its mean satisfies

`f_(n,c_N)(q-y0) >= c/sqrt(n)`.

The density changes only by `o(1)` across the interval of length
`M=o(sqrt(n))`.  Consequently

`K_(c_N)^n E_T(q) >= c M/sqrt(N)`.                                (1)

At the endpoint used in the project,
`f_N(q) asymp log N/N`.  Thus a hoped-for conclusion based only on
`||E_T||_1=O(log T)`, namely

`T^(N-T)E_T(q)=O(f_N(q) log T)`,

is false: the ratio of (1) to the right-hand scale is at least a constant
multiple of `sqrt(N)/log N` (for `T=N^gamma`, the two `log T` factors
cancel).  This remains a counterexample with a height cap and with support
in the late natural tube; neither of those extra properties repairs it.

The mechanism is not a reference-bridge issue.  The later controller
chooses a mean loss `(d_N-T delta)/n=Theta(sqrt(log N/N))`, thereby moving
the required lower-tail terminal point to the *center* of its own sum.

## What a genuine transfer lemma must contain

An unweighted area estimate cannot initialize the r037 late barrier.  A
sufficient replacement must retain a quantitative price for this deliberate
future mean loss.  One possible precise interface is a killed/entropy
estimate: for every future predictable policy `pi`, with
`D_pi=sum_{j>T}(mu-m(c_j))`, one needs a bound of the schematic form

`K_pi^(N-T) E_T(q) <= C f_N(q) log T + R_T(pi)`,

where `R_T(pi)` is absorbed by the late barrier and includes a penalty at
least of reference likelihood size when `D_pi` is near `d_N`.  Equivalently,
one may prove a comparison of the *initial error cone* with a Doob-killed
reference measure, not merely its Lebesgue mass.  A fixed action with
`D_pi=d_N+O(T delta)` has central density `N^(-1/2)`, but its likelihood
relative to the matched reference at that endpoint has the required
moderate-deviation cost of order `exp(-d_N^2/(2 mu N))`.

Concrete next lemma to target: construct a nonnegative error majorant
`E_T<=H_T` together with a predictable compensator `A_pi` such that, for
all policies,

`E_pi[H_T(q-S_(N-T)) exp(-A_pi)] <= C f_N(q) log T`,

and prove separately that the barrier residual pays
`1-exp(-A_pi)`.  Without this paired compensator (or an equivalent
pointwise time-shift/Doob-cone domination), the counterexample above rules
out an initialization based solely on area.

No sharp endpoint assertion follows from this obstruction.
