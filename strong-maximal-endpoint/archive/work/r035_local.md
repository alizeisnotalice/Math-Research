# r035: matched-reference local all-action gain module

Date: 2026-11-06.  This records the large-time local derivation.  Times
below a fixed `t0(A)>e` belong to the short-time initialization, not to the
Fourier tube statement.

Let `a=a_N` solve `2a-m(a)=q/N`, put

`b=m(a)-a>0`, `ell=2a-m(a)=m(a)-2b`,

and let `v_t^a=K_a^t phi_N`, with density `f_t` and the exact derivative
quantity

`d_t(r)=R_t(r)+Q_t(r)/N=-partial_r v_t^a(r)`.

At the critical matching scale, `b=O(sqrt(log N/N))` and
`mu-m(a)=O(b^2)`.  Put `x=r-tm(a)` and choose the common signed tilt `eta`
with tilted mean `r/t`.  On

`|x|<=A sqrt(t log t)`, `t0(A)<=t<=N`,

one has `eta=x/(v_a t)+O(x^2/t^2)`, `|eta|<=C sqrt(log t/t)`,
`d_t>=c f_t`, and the global tilted derivative bounds needed for Taylor
remainders.

## Exact parameter expansion

For

`A_t=partial_c K_c v_{t-1}^a(r)|_{c=a}`,
`B_t=partial_c^2 K_c v_{t-1}^a(r)|_{c=a}`,

the useful form is

`A_t=p_a[ f_t(r)(r-t ell)/(2t) + b Q_t(r)/N + E_A]`,

`B_t=-p_a d_t(r)+E_B`,                                           (A)

with

`|E_A|<=C f_t(r)[(|b|+|eta|)^2+t^-1]`,
`|E_B|<=C[(|b|+|eta|)d_t(r)+f_t(r)/t]`.                           (B)

Derivation: the exact identity is
`A_t=p_a[K_a v_{t-1}^a(r)-v_{t-1}^a(r-a)]` and
`B_t=p_a(2p_a-1)[K_a v-v(r-a)]+p_a v'(r-a)`.
Taylor-expand the first difference around `r-m(a)`.  Its linear part is
`b d_t`; its variance part is `v_a eta f_t/2`; use
`v_a eta=x/t+O(eta^2)` and
`d_t=f_t+Q_t/N+O((|eta|+t^-1)f_t)`.  Since
`b f_t+x f_t/(2t)=f_t(r-t ell)/(2t)`, this gives (A)--(B).
The global shifted density derivative estimate makes the Taylor remainders
integrable against the exponential action tail.  No assertion that every
`b` contribution has a factor `t^-1` is used.

## Global action localization and gain

The r033 action residual is strictly decreasing in `c`; hence its unique
global root is the only maximizer.  Equations (A)--(B) imply

`|c_*(r)-a|<=C(|b|+|eta|+t^-1)`

on the tube, so the compact-action third-derivative estimate applies to
the actual global maximizer.  Taylor's theorem and Young's inequality give,
for every fixed `eps>0`,

`sup_c[K_c v_{t-1}^a(r)-K_a v_{t-1}^a(r)]
 <=(1+eps)(p_a/8) f_t(r)[(r-t ell)/t]^2
   +C_eps f_t(r)[(|b|+|eta|)^2+t^-1]^2
   +C_eps b^2/N`.                                                (G)

The separate `b Q_t/N` in (A) is essential.  Its square divided by the
exact denominator `d_t` is bounded by `b^2 Q_t/N<=b^2/N`; treating
`d_t` as `f_t` too early loses this bound on the lower side of the tube.

## Reference summability

On the natural tube,

`eta^2<=C log(t)/t`,

and the matched critical parameter obeys

`b^2<=C log(N)/N<=C log(t)/t` for `e<=t<=N`.

The last inequality is `t log(N)/N<=log(t)`.  Hence the middle error in
(G) is pointwise at most

`C f_t(r) (log(t))^2/t^2`.

Extend this nonnegative majorant globally.  Exact reference semigroup
propagation gives `P^(N-t)f_t=f_N`, so its full sum is simply

`C f_N(q) sum_{t>=t0}(log(t))^2/t^2=O(f_N(q))`.

No conditional-bridge fourth-moment estimate is needed.  The final
constant error still sums to `O(b^2)=O(log N/N)=O(f_N(q))`, since
`f_N(q) asymp log N/N` at the critical terminal.

Thus (G) gives the desired matched local source plus a reference-summable
error.  It is still only a reference estimate; proving the analogous
all-control MAX correction remains the separate r033/r034 gap.
