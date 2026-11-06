# r032: a sufficient tail-to-endpoint lemma for arbitrary profiles

Date: 2026-11-06.  This note isolates what an all-control CDF theorem must
provide in order to turn a feedback displacement into an **upper** bound for
the feedback endpoint constant.  It does not construct that CDF theorem.

Let `X=X_N^pi` be the terminal sum under an arbitrary predictable policy,
let

`F_pi(s)=P_pi(X<=s)`,
`R_pi(s)=E_pi[e^(X-s) 1_{X<s}]`,

and let `Phi=Phi_N`.  For a nondecreasing `h>=0` set
`J(h)=int_0^infty e^-t Phi(h(t))dt`.

## Profile budget lemma

Put `B_N=sqrt(N) log N`.

If `J(h)<=e^-s` and `E_pi h(X)>1`, then

`F_pi(s) <= R_pi(s)
 + N^-1 E_pi(X-s)_+
 + C N^-2 E_pi(X-s)_+^2
 + Theta_pi,h,N(s)`,                                             (P)

where

`Theta_pi,h,N(s)=E_pi[(h(X)-1)_+ 1_{X>s+B_N}]`.

Proof: monotonicity gives `e^-t Phi(h(t))<=J(h)`, so below `s`,
`h(t)<=e^(t-s)`.  On `0<=u<=B_N`, inversion of
`y+(N-1)log(1+y)<=u`, `y=log h`, gives
`(h(s+u)-1)_+<=u/N+C u^2/N^2`.  Above the band, the same pointwise inverse
is not used: the far-tail term is deliberately retained with the actual
profile `h`.  Finally

`1<Eh(X)<=R_pi(s)+1-F_pi(s)+E[(h(X)-1)_+1_{X>=s}]`,

which yields (P).  This is valid for every layered/soft/capped profile, not
only the plateau used for a lower bound.

## A single, nonvacuous endpoint interface

Define the piecewise terminal payoff

`H_{N,s}(x)= exp(x-s)` for `x<s`,
`H_{N,s}(x)=1+(x-s)/N` for `x>=s`.

For every policy, conditional means satisfy `E[Z_c|past]<=mu` and
conditional variances are uniformly bounded.  Writing
`X=N mu-D+M` with `D>=0` predictable drift loss and `M` a martingale gives,
for `s=N mu-d_N+O(log N)`,

`sup_pi E_pi(X-s)_+^2=O(N log N)`.                              (M)

The near-band calculation in (P), the tail lemma below, and (M) yield the
uniform **arbitrary-profile reduction**

`J(h)<=e^-s  =>  E_pi h(X)
 <=E_pi H_{N,s}(X)+O(log N/N)`,                                  (R)

for every nondecreasing `h` and every predictable policy.  This reduction
keeps the drift effect exactly; indeed

`E H_{N,s}(X)=1-F_pi(s)+R_pi(s)+N^-1E_pi(X-s)_+`
`=1-F_pi(s)+R_pi(s)+(d-D_pi+L_pi(s))/N`.

Thus neither a false all-policy `L_pi=O(1)` bound nor a separate
all-policy density comparison is required.

Let `s_Q=N mu-d_N+a_Q log N+o(log N)`,
`a_Q=(9-mu)/24`, and `g=(mu-1)/8`.  The exact remaining sufficient theorem
is the following one-payoff all-control estimate: for every fixed
`epsilon>0`, at

`s=s_Q+(g+epsilon)log N`,

`sup_{predictable pi} E_pi H_{N,s}(X)
 <=1-c_epsilon f_ref,N log N+o((log N)^2/N)`,                   (E)

where `c_epsilon>0` and `f_ref,N asymp log N/N`.  Equivalently the
right-hand deficit is a positive constant times `(log N)^2/N`.

Then (R) contradicts `E_pi h(X)>1`, because its reduction error is only
`O(log N/N)=o(f_ref,N log N)`.  Therefore

`log R_N<=s_Q+g log N+o(log N)
 =N mu-d_N+((3+mu)/12)log N+o(log N)`.

This formulation says exactly what an all-control Poisson/CDF argument must
deliver.  A CDF estimate alone is useful only if it is strong enough to
imply (E) for this exponential/linear payoff; an upper CDF bound for one
constructed policy still proves only a lower bound for `R_N`.

## Uniform Fenchel tail absorption (actual budget, not pointwise inverse)

This part is available uniformly for all policies.  The conditional one-step
density obeys `q_c(z)<=2e^-z 1_{z>0}`.  Iterating conditional densities,
including predictable choices of `c`, gives

`f_pi,N(t)<=2^N e^-t t^(N-1)/Gamma(N)`.                         (F)

A pointwise inverse majorant is false even for the reference Gamma law:
with fixed `N`, `Phi_N^(-1)(e^(t-s))` grows like
`e^(t-s)/(t-s)^(N-1)`, while its product with
`e^-t t^(N-1)/Gamma(N)` tends to a nonzero constant.  Its tail integral
therefore diverges.  The following Fenchel step is essential.

Fix `s=N mu+o(N)`, take `epsilon_N=exp(-kappa(log N)^2)`, and set the
Fenchel multiplier `Lambda=epsilon_N e^s`.  Young--Fenchel on the tail
gives

`int_{t>s+B_N}h(t)f_pi,N(t)dt
 <=Lambda J(h)
   +Lambda int_{t>s+B_N}e^-t
      Phi_N^*( e^t f_pi,N(t)/Lambda )dt`.                       (Y)

The first term is at most `epsilon_N`.  Let `h_*(q)` maximize the conjugate.
From (F), Stirling, `s=N mu+o(N)`, and `-log epsilon_N=o(N)`,

`q(t)^(1/(N-1))
 =(e^t f_pi,N(t)/Lambda)^(1/(N-1)) <=C t/N`.

For `q>N`, the identity for `Phi_N'` implies
`h_*(q)<=exp(Ct/N)`; for `q<=N` it is at most one.  Since
`Phi_N^*(q)<=q h_*(q)`, the second term in (Y) is at most

`E_pi[exp(CX/N)1_{X>s+B_N}]`.

Uniformly for small positive `theta`, the conditional mgf satisfies
`E[exp(theta Z_c)|past]<=exp(mu theta+C theta^2)`.  Choose
`theta` of order `(log N)/sqrt(N)` in the Chernoff estimate, absorbing the
fixed `C/N` weight.  Because
`s+B_N-N mu=sqrt(N)log N-d_N+O(log N)` is positive of order
`sqrt(N)log N`, this last expectation is `<=exp(-c(log N)^2)`.

Thus

`sup_{pi, J(h)<=e^-s} E_pi[h(X)1_{X>s+B_N}]
 <=exp(-c'(log N)^2)`,                                         (T)

and hence the same holds for `Theta_pi,h,N(s)`.  This proves the tail term
used in (R).  Unlike the invalid pointwise inverse majorant, (T) uses the
single integrated Orlicz budget `J(h)`.

## What is essential, and what is not

The deficit scale in (E) is natural: a shift of fixed multiples of `log N`
changes the critical reference CDF by `f_ref,N log N asymp(log N)^2/N`.
Errors merely `o(1)` or even `O(sqrt(log N/N))` are insufficient.  The
Fenchel argument closes the high-layer issue for arbitrary `h`; the genuine
remaining task is the all-control estimate (E) for the one explicit payoff.

## Formal r032 report audit (sections 2--4)

The formal report states the scopes correctly.  Its Fenchel estimate is
uniform over every predictable policy and every nonnegative measurable
budget input; it deliberately does not use monotonicity.  The elementary
density bound survives feedback by conditional simplex integration, and the
positive-mgf Chernoff step has a uniform action range including `c=infty`.
The conjugate argument is valid: for `a>N`, the maximizer has `h=e^y>1`
and `Phi_N'(e^y)>=y^(N-1)`, while (2), Stirling, and the chosen multiplier
give `y<=Ct/N`.

Monotonicity enters only in the separate reduction to `H_{N,s}`.  The
martingale decomposition has predictable nonnegative drift loss and
uniform conditional variance, yielding the stated `O(N log N)` overshoot
bound.  Formula (9) retains, rather than assumes away, both drift and
low-side terms.  Finally, (10) is explicitly a pending all-control Bellman
estimate, so the conditional consequence (11) is not reported as a theorem.
