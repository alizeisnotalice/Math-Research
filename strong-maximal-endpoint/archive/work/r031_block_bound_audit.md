# r031: 2/3-block KL bound audit — 2026-11-06

## Verdict

The proposed bounded-block route is valid for `c,a` in a compact interval
`J=[a0,a1]` with `a0>0`, provided the segment between `a,c` stays in `J`.
It gives a uniform `D_b(c||a)<=C_J b(c-a)^2` for every `b>=2` (with the
usual 2/3 decomposition), and the fixed-`b` Fisher coefficient is
`I_b(c)/2`.  No large-`b` Fisher asymptotic is used.

## Scores and knots

Put `p_c=(1+e^c)^-1`.  The derivative identity

`partial_c f_{b,c}(t)=b p_c [f_{b,c}(t)-f_{b-1,c}(t-c)]`

gives `score_b=b p_c(1-f_{b-1,c}(t-c)/f_{b,c}(t))` a.e.  For `b=2`,

`P_2=t+2(t-c)_+ +(t-2c)_+`,

and, when the shifted term is nonzero,

`f_{1,c}(t-c)/f_{2,c}(t) <=(1+e^c)/c`.

Thus

`|score_2|<=M_2(c):=max(2p_c,2/c-2p_c)`.

For `b=3`, `P_3>=t^2` and `P_2(t-c)<=4t`; hence

`f_{2,c}(t-c)/f_{3,c}(t)<=8(1+e^c)/c`,

which gives

`|score_3|<=M_3(c):=3p_c+24/c`.

These bounds include the regions below the first knot, where the shifted
term is zero.  At `t=kc` the classical derivative may have one-sided
values, but the density is positive for `t>0` because of the unshifted
Gamma component and the finite knot set has zero Lebesgue measure.  It is
therefore harmless for score/KL integrals.

## Likelihood-ratio and KL bound

Let `M_J=max_{c in J}(M_2(c) vee M_3(c))`.  Absolute continuity in the
parameter for a.e. `t` gives

`|log(f_{r,c}/f_{r,a})|<=M_J|c-a|`, `r=2,3`.

For `x=log(f_{r,c}/f_{r,a})`, use

`e^x x-e^x+1 <= (1/2)e^{|x|}x^2`

and `E_a(e^x-1)=0`.  This proves exactly

`D(f_{r,c}||f_{r,a}) <= (1/2)e^{M_J|c-a|}M_J^2(c-a)^2`,

so the asserted quadratic upper bound does not require a hidden Fisher
calculation.

For `b=2m` write the sum as `m` independent 2-blocks; for
`b=2m+3` use `m` 2-blocks plus a 3-block.  Data processing under the sum
map gives `D_b<=mD_2(+D_3)<=C_Jb(c-a)^2` (the statement should explicitly
say `b>=2`; `b=1` needs its own elementary treatment if ever used).

## DQM / fixed-b expansion

The uniform score bounds give `partial_c sqrt(f)=score*sqrt(f)/2` in
`L^2`.  Away from finitely many moving knots the pointwise difference
quotient converges; the bounds above dominate it.  The moving knot strips
have measure `O(|delta|)` and bounded difference quotients, hence vanish
in the `L^2` limit.  Thus the family is differentiable in quadratic mean
for fixed `b=2,3`, and by finite convolution for every fixed `b`.  Therefore

`D(f_{b,c+delta}||f_{b,c})=(I_b(c)/2)delta^2+o(delta^2)`,

where `I_b(c)=int score_{b,c}^2 f_{b,c}`.

## Local computational check

MATLAB evaluated the explicit scores on a dense `t` grid for
`c=0.2,0.7,mu`; the ratios `max|score_2|/M_2` were
`0.998996,0.999548,0.999757`, and `max|score_3|/M_3` was below `0.034` in
all three cases.  This supports but does not replace the algebra above.

## Joint-limit obstruction (also accepted)

The fixed-`b` Fisher expansion cannot be made into a `b`-uniform
`D_b<=C(c-mu)^2` statement.  The entropy variational inequality gives

`D(f_{b,c}||f_{b,mu}) >= b sup_lambda
 {lambda m(c)-log M_mu(lambda)}`.

Since `m'(mu)=0` and `m''(mu)=-p_mu`, for `c=mu+h` one has
`m(c)=mu-p_mu h^2/2+O(h^3)`.  Taking
`lambda=-p_mu h^2/(2mu)` gives

`D_b >= b[p_mu^2 h^4/(8mu)+O(h^5)]`.

At `h=b^(-1/4)` this has a positive limiting lower bound while `h^2` tends
to zero.  It is fully compatible with the earlier `C_J b h^2` block upper
bound; it only says that the Fisher-radius must shrink with `b` in any
joint-limit claim.  MATLAB finite differences confirmed
`(m(mu+h)-mu)/h^2 -> -p_mu/2` and gave
`p_mu^2/(8mu)=0.0046385661825`.

## Final review of r031 sections 3--5

**Accepted.**  In the common tilt identity the numerator is centred at
`d=K_mu'(eta)-mu`, while the denominator is at zero.  The exact untilting
factor is `exp(eta mu)/M_mu(eta)`, so there is no shift or sign loss in
(12).  The `A_eta/k` Edgeworth term is the same single-increment term for
`b` and `b-1`; its cancellation in their ratio is valid.  The odd first
Gaussian correction vanishes at zero, and the displayed global second
derivative bound gives (15), hence the required relative `1/b` ratio.

On the `sqrt(log b)` central band, (17) has an error whose product with the
polynomial score remains uniformly negligible.  Off the band, the earlier
global `|score_b|<=C_J b` bound combined with a Chernoff tail of order
`b^-M`, `M>3`, makes the Fisher tail `o(1)`.  The exact fourth moment of
the standardized sum gives the limit 2, and the uniform eighth-moment bound
legitimately removes the small-probability tail by Cauchy--Schwarz.  Thus
the proof does not integrate a merely pointwise CLT, and it proves
`I_b(mu)->p^2/2`.
