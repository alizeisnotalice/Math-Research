# r034: annular reference sources stay `O(f_N)` at the critical terminal

Date: 2026-11-06.  This is a reference-convolution estimate only.  It does
not bound the all-control Bellman correction.

Let `S_t` be the independent common-`mu` reference sum with density `f_t`.
Fix a critical terminal window

`q=N mu-d_N+O(log N)`,

where `d_N/sqrt(N log N)->sigma`.  Put `L=log N` and assume sources are
active only for `t0<=t<=N`, where `t0>=C0 L`.  Suppose, for fixed exponents
`P,Q`, that

`|e_t(r)|<=C N^P t^Q f_t(r)
  1_{|r-t mu|>=A sqrt(t L)}`.                                  (S)

The same result applies to a sum of finitely many such sources.

## Conditional bridge bound

Use the common signed tilt `eta` satisfying `m_eta=q/N`.  It is uniformly
`O(sqrt(L/N))`; after tilt all exponential rates and all fixed moments are
uniform.  Exact cancellation of tilt factors gives, for `t<N`,

`P(S_t in dr | S_N=q)
 =f_t^eta(r) f_{N-t}^eta(q-r)/f_N^eta(q) dr`.                    (B)

At its mean the denominator is at least `c/sqrt(N)`.  Every nonempty
suffix density has a uniform `L^infty` bound (use the exponential factor),
including `N-t=1`.  Hence

`P(E | S_N=q)<=C sqrt(N) P_eta(S_t in E)`.                       (C)

This deliberately remains valid when `t` is within one step of `N`.

For large `N`, uniformly in the stated terminal window,

`sqrt(t)|m_eta-mu| <=(sigma+delta)sqrt(L)`

for any fixed `delta>0`.  If `A>sigma+2delta`, the annulus event in (S)
forces a centred tilted deviation at least

`(A-sigma-2delta)sqrt(t L)`.

Since `t>=C0 L`, uniform Bernstein/Chernoff for the tilted increments gives

`P_eta(|S_t-t mu|>=A sqrt(tL))
 <=2 exp[-c0(A-sigma-2delta)^2 L]`.                              (D)

For `t=N` the event is empty for all large `N`, because
`|q-Nmu|<A sqrt(NL)`.  Thus (C)--(D) cover all active times, including the
near-terminal range.

## Propagation estimate and threshold for `A`

The exact bridge identity for a source `e_t=f_t rho_t` is

`(P^(N-t)e_t)(q)
 =f_N(q) E[rho_t(S_t)|S_N=q]`.                                  (E)

Combining (S)--(D), then summing over `t`, yields

`sum_{t=t0}^N |(P^(N-t)e_t)(q)|
 <=C f_N(q) N^(P+1/2+beta_Q-c0(A-sigma-2delta)^2)`,               (F)

where `beta_Q=max(Q+1,0)` (at `Q=-1`, insert an inessential `log N`
factor).  Therefore it is sufficient to choose

`c0(A-sigma-2delta)^2 >= P+1/2+beta_Q`

to obtain an `O(f_N(q))` terminal error, and strict inequality gives
`o(f_N(q))`.  This states explicitly how the annulus width must dominate
the polynomial prefactor.  If the source is `C t^-K f_t` on the annulus,
take `Q=-K`; the exact harmonic/power sum improves (F) in the evident way.

## Scope and limitations

The estimate needs hard annular support (or a separately verified tail
decomposition).  A bare fixed Gaussian factor
`exp(-a(r-tmu)^2/t)` is not automatically interchangeable with hard support:
its bridge expectation must be evaluated, and can retain central mass.  The
small-time regime `t<t0` is likewise not covered by (D); it belongs to a
short-time initialization/splicing argument.  No all-control conclusion or
sharp endpoint coefficient follows from this reference estimate.
