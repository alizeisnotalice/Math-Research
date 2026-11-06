# r016 independent theory audit — 2026-11-06

Scope: only `rounds/r016/report.md`; code and Git untouched.

## Verdict

Accepted. The stated uniform, fixed-window CDF/quantile theorem follows from
the displayed arguments. It is correctly scoped to independent fixed-
proportion inputs and does not claim a feedback or endpoint-constant result.

## Checks of the critical links

1. For
\[
M_c(\lambda)=\frac{1+e^{-(1+\lambda)c}}
 {(1+\lambda)(1+e^{-c})},\qquad p_c=(1+e^c)^{-1},
\]
direct differentiation gives exactly
`d_c M_c=p_c(M_c-e^{-lambda c})`. Laplace inversion therefore yields
\[
 \partial_cF_{n,c}(q)=np_c\{F_{n,c}(q)-F_{n-1,c}(q-c)\},
\]
and the sign in the implicit derivative (3) is correct. A local Wolfram 15
symbolic check returned zero for this one-step identity.

2. The Fourier lemma is sufficient and genuinely self-contained. The
exponential component gives the global integrable bound (4); after extracting
half the powers it gives uniform exponential control outside the small-
frequency interval. On that interval the uniform cubic log-characteristic
remainder and the negative real quadratic part give (6). Its weighted
integrals have orders `N^-1` (density at zero) and `N^-3/2` (first derivative,
after the odd Gaussian integral is cancelled). This validates precisely the
otherwise delicate assertion `g_N'(0)=O(N^-3/2)`.

3. The tilt is feasible because `a=q_c-c` is eventually in `(0,Nm(c))`;
the FGL localisation only needs `x_c=sqrt(log n)(1+o(1))`, uniformly in a
bounded `t` window. Thus `lambda=rho/sigma(1+o(1))` is adequate for the
claimed leading coefficient.

4. The Taylor bookkeeping is correctly centred at `X=Z_c-c`, not at an
absolute moment. In particular,
\[
 E X=m(c)-c=O(\rho),\quad
 f'_{N,c}(a)/f_{N,c}(a)=\lambda+O(N^{-1}),
\]
and the exponential-weighted third moment is uniformly bounded. Hence (11)
has the stated remainder. Repeating this for the density retains the first
moment, yielding
`f_{n,c}(q_c)/f_{N,c}(a)=1+O(rho*lambda+lambda^2+1/N)=1+O(L/n)`.
There is no lost order here. The extension at zero has the needed `C^2`
regularity for large `N`: each shifted Gamma component has shape `N` and
vanishes to sufficiently high order at its left endpoint.

5. Substitution in (3) gives
\[
 q'_c/(n\rho)=-p_*\,(t+\sigma/2)+o(1)
\]
uniformly for bounded `t`. Since `dc=rho dt` and `n rho^2=log n`, integration
does prove (1), including uniformity. The maximum is at `t=-sigma/2` and has
value `p_*sigma^2/8=(mu-1)/8`, so (2) follows.

Actual local calculation: Wolfram 15 returned zero for the Laplace derivative
identity and evaluates `(mu-1)/8` as
`0.0348080678451342243887`, matching the stated constant.

No blocker found in §3's outer-frequency integration or §5's remainders and
density-ratio argument.
