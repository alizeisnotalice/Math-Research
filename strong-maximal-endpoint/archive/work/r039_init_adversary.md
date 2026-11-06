# r039: finite-time global initialization audit

Date: 2026-11-06.  The proposed statement is credible for each fixed
`k>=4` (and all sufficiently large matched `N`), but its right-tail proof
must use an exponential-moment screen rather than the uniform bound on the
payoff.

Let

`I_k=v_k^a exp((alpha J_k+s0 f_k)/v_k)+Lambda k b^2/N`,

where `J_k` is the nonnegative quadratic Poisson potential.  Since
`exp(x)>=1+x`, it suffices in the central/small regions to dominate the
true-minus-reference difference linearly.  The three regions are as follows.

## Left and small positive states

For `r<=0`, every increment is nonnegative and the terminal payoff is
affine.  Hence exactly

`U_(k,N)(r)-v_(k,N)^a(r)=k[mu-m(a)]/N=k delta/N`.

Choose fixed `Lambda` with `Lambda b^2>=delta`; the floor dominates this
identity.  For `0<r<=r0`, write

`phi_N(r-S)=1+(S-r)/N-psi_N(r-S)1_{S<r}`,

with `0<=psi_N(u)<=u`. Dropping the penalty for the MAX value gives

`U_k-v_k^a<=k delta/N+C r P_a(S_k<r)<=k delta/N+C r^(k+1)`.

Uniformly for matched `a` in a compact neighbourhood of `mu`, the all-zero
Bernoulli component gives `f_k^a(r)>=c r^(k-1)` on this interval.  A fixed
`s0` therefore absorbs the `r^(k+1)` term, and the same floor absorbs the
drift term.  This also rules out an `r=0` counterexample.

## Compact middle region

After fixing `0<r0<R`, the family of matched actions is compact and
`inf_{a,r in [r0,R]} f_k^a(r)>0`.  The true payoff is uniformly bounded
(`phi_N(r-S)<=1+S/N` and every action has mean at most `mu`), so one fixed
additional amount in `s0 f_k` dominates `U_k-v_k^a` on this compact region.
This is a legitimate existential, not quantitative, choice of `s0`.

## Correct far-right argument

The naive comparison with `U_k<=1+k mu/N` fails because the proposed
`I_k` can still tend to zero when the linear coefficient in `J_k/v_k` is
below one.  Instead, for every `theta in [1/N,1)`,

`phi_N(u)<=exp(-theta u)` for all real `u`,

so

`U_(k,N)(r)<=M_*(theta)^k exp(-theta r)`.                           (1)

The finite-step constants are uniform in matched `a`: for large r,

`f_k^a(r)>=c e^(-r)r^(k-1)`,
`v_k^a(r)<=C e^(-r)r^k`,
`J_k(r)>=c f_k^a(r)(r-k ell)^2>=c e^(-r)r^(k+1)`.

The last inequality uses only the last nonnegative source in the Poisson
recursion. Therefore `alpha J_k/v_k>=kappa r` for a fixed
`0<kappa<1`, after increasing `R` if necessary. Hence

`I_k(r)>=c e^(-(1-kappa)r)r^k`.

Choose `theta=1-kappa/2`; for large `N` it belongs to `[1/N,1)`, and (1)
is exponentially smaller than this lower bound.  Thus the far right is
also uniformly dominated.  Wolfram verified the Bernoulli exponential
coefficient `E exp(c B)=2/(1+e^(-c))`, used in the corresponding tail
constants.

## Exact remaining task / next useful boundary

The proof above supplies a fixed-time initialization majorant, but not the
late `N^gamma` initialization required by r037.  A promising next
fixed-time one-step boundary is a moderate-deviation comparison with
`lambda=sqrt(log t/t)` on the *negative* side: the small-state inequality
and the r036 local MAX lemma meet there.  What must be shown is a uniform
overlap estimate for the actual finite-`k` majorant and the late barrier;
no global fixed spatial shift can supply it because r038's linear shift
obstruction applies.  This is not yet proved here.
