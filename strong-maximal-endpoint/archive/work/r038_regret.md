# r038: exact regret-killed transfer interface

Date: 2026-11-06.  This is a feedback-respecting replacement for an
unjustified area-to-reference-density propagation claim.

Let the late r037 barrier be `W_t` and define its exact one-step regret

`R_t^c(r)=W_t(r)-K_cW_(t-1)(r) >= 0`.                              (1)

For a predictable backward policy, write `r_t=q-sum_(j=t+1)^N Z_j`.
Conditional expectation gives the exact telescope

`W_N(q)-E_pi W_T(r_T)=sum_(t=T+1)^N E_pi R_t^(c_t)(r_t)`.           (2)

There is no reference-path replacement in (2).

## A proved remote-action killing bound

On a left-exterior state where the selected coarse branch is specifically
`W_t=O_t^-`, use `W_(t-1)<=O_(t-1)^-` and r034's affine/Gamma barrier:

`R_t^c(r) >= O_t^-(r)-K_cO_(t-1)^-(r)
            >= [mu-m(c)]/N`.                                      (3)

The last inequality follows because the Gamma correction is nondecreasing
and `Z_c>=Exp(1)` under the common coupling.  A local Wolfram algebra check
gave the affine difference exactly as `(mu-m(c))/N`.

Consequently every policy satisfies

`sum E R_t >= (1/N) E sum [mu-m(c_t)]
                    1_{r_t in left exterior, W_t=O_t^-}`.          (4)

This is a genuine mean-loss absorption statement.  It is deliberately
restricted to the selected left coarse branch; replacing that condition by
the unproved assertion `W=O^-` everywhere would be invalid.

## Minimal killed-error envelope: a precise next lemma

Let `E_T=(U_true,T-W_T)_+`.  Define recursively, with all feedback retained,

`e_T=E_T`,
`e_t(r)=sup_c [K_c e_(t-1)(r)-R_t^c(r)]_+`,  `T<t<=N`.               (5)

This is the smallest nonnegative sequence satisfying

`K_c e_(t-1)(r) <= e_t(r)+R_t^c(r)` for every action.               (6)

Iterating (6) along an arbitrary predictable policy and combining with
(2) gives

`E_pi E_T(r_T) <= e_N(q)+sum E_pi R_t`,
`E_pi[W_T(r_T)+E_T(r_T)] <= W_N(q)+e_N(q)`.                          (7)

Taking the supremum over policies shows the exact sufficient endpoint
condition

`e_N(q) <= C f_N(q) log T`.                                        (8)

Under (8), the late barrier plus the true initialization is transferred
with the desired error.  Thus (8), for the *actual structured* `E_T`, is a
well-posed next lemma.  It is neither an L1 estimate nor a fixed-reference
bridge assertion; its recursion includes the same feedback and the same
regret that produced the late barrier.

## Limitation

The bound (4) alone cannot prove (8) from an area estimate: a future policy
can distribute a total mean loss of order `d_N` while centering the terminal
convolution.  Its affine killing is only `d_N/N`.  Therefore a proof of
(8) must use further pointwise structure of `E_T` (or a stronger Doob/entropy
killing term in the interior), not just its mass.  This is a scope warning,
not a new endpoint claim.

## Strengthening: the left coarse branch is selected on the whole left exterior

For `r<=tm(a)-A sqrt(t logt)`, put `x=t mu-r`. Then
`x>=A sqrt(t logt)`. Jensen gives
`M_*(theta)>=M_mu(theta)>=exp(mu theta)`, so the positive exponential
screen obeys

`O_t^+(r)>=inf_(theta>=1/N) exp(theta x)=exp(x/N)`.                 (9)

The affine part of `O_t^-` is `1+x/N`. If `r<=rho_N`, its Gamma correction
vanishes. If `r>rho_N`, then, since `t<=N`,

`r-rho_N<=mu t-N+1<.3t`

for all sufficiently large `N`. On `r<=mu t`, `b_N` is uniformly bounded,
and the Gamma convolution formula consequently gives

`0<=B_t^Gamma(r)<=C P(Gamma(t,1)<.3t)<=C exp(-ct)`.                (10)

On the other hand

`exp(x/N)-1-x/N>=x^2/(2N^2)>=c t logt/N^2`.

For `t>=N^gamma` the last quantity dominates `exp(-ct)` uniformly. Hence
`O_t^-(r)<=O_t^+(r)` throughout the whole late left exterior. Combined
with r037's `I_t>=O_t^-` there, this proves exactly `W_t=O_t^-` on that
region, rather than merely that some coarse minimum is selected. Therefore
the killing estimate (3)--(4) applies to every late left-exterior state.
The argument has no consequence for the terminal initialization problem by
itself.
