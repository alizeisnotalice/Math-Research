# r059: shifted response tails and the natural-angle continuation

## Scope and conclusion

This note concerns the frozen symmetric cap-two source of r057/r058, with
`N=2b`, `L=log N`, a fixed angular radius `A`, and the regular KKT hypotheses
already used there.  It fills the *axis issue* which prevents copying the
one-block density-supremum proof from r058 after a receiver translation.  It
also records the resulting continuation argument.  It does not prove a sign
for `K`, a strict endpoint gain, or any assertion about the general strong
maximal constant.

The result is conditional only on the already-established r054 compact tilt
bridge, the r057/r058 finite list of score/residual bounds, and the standard
regular differentiation/coarea interchanges.  No global bound on a shifted
one-block pushforward density is asserted (indeed it is false in general).

Put `u=(d+v,d-v)` and write, in source coordinates,

\[
 \widetilde f_i(t)=f_{b,a}(t-u_i){\bf1}_{t>u_i},\qquad
 \mathcal R_i(t)=\frac{\widetilde f_i(t)}{f_{b,c}(t)}.
\]

Thus a shifted response form is the source complete response form with its
`f_{b,c}(t_1)f_{b,c}(t_2)` measure multiplied by
`\mathcal R_1\mathcal R_2`; call it `B_R`.  This convention includes the
lower coarea term and the soft term.  The exact derivatives, with source
scores `k_i=partial_t log R_{b,c}` and target common scores
`\widetilde s_i=partial_a log R_{b,a}(t-u_i)`, are

\[
 G_d=B_R(1,k_1+k_2),\quad G_v=B_R(1,k_1-k_2),
\]
\[
 G_{ad}=B_R(k_1+k_2,\widetilde s_1+\widetilde s_2),\quad
 G_{av}=B_R(k_1-k_2,\widetilde s_1+\widetilde s_2).                 \tag{1}
\]

Here and below identities include the lower jump distributionally.  The
formula is obtained by changing the receiver integral to source coordinates;
the derivative of the frozen profile supplies the source `k`, whereas the
receiver-action derivative supplies the target `s`.

## The support estimate that removes the axis singularity

For `q=e^c`, the binomial/Gamma formula has the useful cancellation

\[
 R_{b,c}(t)=\frac{(q/(1+q))^b}{\Gamma(b)}
 \sum_{j=0}^b {b\choose j}(t-cj)_+^{b-1}.
\]

Consequently, uniformly for the compact `c` window,

\[
 R_{b,c}(t)\le {2^b t^{b-1}\over\Gamma(b)}.                       \tag{2}
\]

The response is supported on the lower level `lambda R_1R_2=1` or on the
soft strip `N<lambda R_1R_2<Phi_N'(2)`.  Since `H=1` is feasible, `P<=1`,
and r054 gives `lambda/P<=2N`, one has `lambda<=2N=4b`.  Hence every point
that can contribute to either part of `B_R` satisfies

\[
 t_1t_2\ge\left\{\frac{\Gamma(b)^2}{4b\,4^b}\right\}^{1/(b-1)}
 \ge {b^2\over128e^2}=:\eta b^2 \qquad(b\ge6).                  \tag{3}
\]

For the last inequality use `(b-1)! >= ((b-1)/e)^(b-1)`,
`(4b)^(-1/(b-1))>=1/2`, `4^(-b/(b-1))>=1/16`, and
`(b-1)^2>=b^2/4`.  This is deliberately a loose constant.  The accompanying
Wolfram receipt checks the binomial cancellation and evaluates the displayed
elementary lower bound; the proof of (3) is the preceding analytic chain.

This is the essential point missing from a naive shifted version of r058:
if `u_i<0`, then `\widetilde f_i(t)/f_{b,c}(t)` can diverge as `t` tends to
zero, but a response point with `t_i` near zero has its other source
coordinate at least `eta b^2/t_i`, far out in a Gamma tail.

## Finite-weight shifted coarea-tail lemma

Fix `A,M`.  There are `D=D(A,M)` and `b_0` such that, if

\[
 |v|\le A\sqrt N,\quad |d|\le C_A\sqrt L,\quad
 |a-c|\le C_A L^{3/2}/N,                                             \tag{4}
\]

then the lower and soft pieces of `B_R` outside

\[
 \mathcal C_D=\{\max_i|t_i-bm(c)|\le D\sqrt{bL}\}
\]

are `O_{A,M}(N^{-M})` for each of the following response weights:

\[
 1,\quad k_i,\quad k_i\widetilde s_j,\quad
 k_i(\widetilde k_j-1),                                               \tag{5}
\]

where `\widetilde k_j=\partial_t\log R_{b,a}(t-u_j)`.  These are precisely
the weights in (1) and in the `d`-derivative used below.  Factors of
`H_0-1` are harmless because they are bounded.  Target residuals and score
derivatives are deliberately excluded from this response lemma: they are
needed only under an ordinary target expectation, for which (9) below is the
appropriate argument.  This is a finite list, not a claim for arbitrary
polynomial weights.

Here is a proof sufficient to expose all formerly singular steps.

1. On `C_D`, the compact local-limit bridge, applied at `c` and `a`, gives

   \[
   \exp(-C_A\sqrt L)\le\mathcal R_1\mathcal R_2
   \le\exp(C_A\sqrt L).                                             \tag{6}
   \]

   The `sqrt L` loss comes from translating a `D sqrt(bL)` local-limit
   coordinate by `O_A(sqrt b)`.  It is only `N^{o(1)}`.  The same bridge gives
   `k_i-1=O(sqrt(L/N))`, `\widetilde s_i=O_A(L)`, and the required bounded
   central versions of the target scores in (5).

2. Parameterize the lower level by `U_i=log R_{b,c}(t_i)`.  Because
   `k_i>=(b-1)/t_i`, its target pushforward is

   \[
    \widetilde q_i(U_i)=\widetilde f_i(t_i)/k_i.
   \]

   It is *not* uniformly bounded: a shifted density ratio can be singular at
   a source axis.  Instead integrate the coarea delta in the other coordinate.
   A typical displayed weight has at most one source score from each
   differentiated frozen profile, and has density bounded by

   \[
   \widetilde f_1(t_1)\widetilde f_2(t_2)
   k_1^{r_1}k_2^{r_2-1},dt_1,\qquad 0\le r_i\le1.                 \tag{7}
   \]

   If `t_1<epsilon b`, choose `epsilon<eta/(4\sup_{c\sim\mu}m(c))`.
   By (3), `t_2>(eta/epsilon)b`; after an `O_A(sqrt b)` shift this is a
   fixed positive large-deviation distance above the target mean.  Each
   negative source power is removed using

   \[
    k_i\le Cb/t_i,\qquad t_1^{-1}\le t_2/(\eta b^2),                \tag{8}
   \]

   where the first bound is the Euler score identity and the second is (3).
   Positive source score factors use the upper bound in (8); reciprocal
   score/Jacobian factors use the lower bound `k_i>=(b-1)/t_i`.
   Thus (7), including all source-score factors, is at most a fixed
   polynomial in `t_2/b` times the target Gamma-mixture density at `t_2-u_2`.
   The source-side density `\widetilde f_1(t_1)dt_1` is retained and
   integrated, rather than discarded: its total mass on the axis interval is
   at most one.  The supremum of the remaining polynomial times the second
   target density on this fixed upper-tail region is `O(e^{-c_A b})`.
   The target scores `\widetilde s_j,\widetilde k_j` are finite sums of fixed
   polynomial and inverse powers of shifted Gamma variables; inverse powers
   lower the shape by a fixed amount.  Thus their weighted Gamma-mixture
   densities preserve this exponential upper-tail bound for large `b`.  The case
   `t_2<epsilon b` is symmetric.

3. On the non-axis set `t_1,t_2>=epsilon b`, the source scores are bounded
   by `C_epsilon`; the target densities times each fixed target score have
   Gamma-shift maxima polynomial in `b` (including the coarea Jacobian
   `1/k_j<=t_j/(b-1)`).  If, say, coordinate one leaves `C_D`, retain its
   target density and integrate its polynomial weight against the target
   law.  Hölder--Chernoff bounds this moment by a polynomial times
   `exp(-cD^2L)`.  Take the supremum only of the other weighted target
   density, which is polynomial in `b`; then choose `D` after this finite
   list.  This covers the unbounded upper region without a global shifted
   density-ratio bound and is uniform in the source lower level.  The other
   coordinate is identical by symmetry.

4. In the soft term the same product support (3) holds and its factor is at
   most `6/N`.  The axis transfer and the non-axis target moments above apply
   unchanged.

This proves the lemma without a fictitious shifted `L^infty` pushforward
bound.  Notice that merely saying “the other coordinate is large” would not
be enough: (8) is what transfers every possible source-axis inverse power to
that coordinate before the Gamma tail is used.

There is a separate, easier companion fact which must not be conflated with
this response lemma.  For the *ordinary* shifted target product measure,
`t_i-u_i` is a Gamma--Bernoulli sum with mean `bm(a)`, while the source
central window has centre `bm(c)` and the displacement between these centres
is only `O_A(sqrt b)`.  Direct Chernoff bounds therefore give

\[
 E_{a,u}[W\,{\bf1}_{\mathcal C_D^c}]=O_{A,M}(N^{-M})                \tag{9}
\]

after enlarging `D`, for every fixed target score/residual weight `W` used
below.  In particular this includes `W=1,|Z_a|,Z_a^2`.  This ordinary-measure
statement does **not** use (3), and it has no source-axis density-ratio
singularity.

## Uniform continuation

First make the bootstrap non-circular.  Fix `B` and work on the maximal
regular interval on which

\[
 |d|\le B\sqrt L,\qquad |a-c|\le B L^{3/2}/N.                     \tag{10}
\]

For fixed `A,D`, the constants in the central local-limit comparison are
independent of `B` once `N` is large: `B sqrt(L)=o(sqrt b)` and
`B L^{3/2}/N=o(1)`.  Thus the constants below depend on `A,D` (and the
original fixed `c` window), but not on this provisional `B`.

Let `(d(v),a(v))` solve `G=1,G_a=0` and start at `(0,c)`.  The response-tail
lemma, (6), and the r057 normalization give throughout (10),

\[
 G_d\ge N^{-1-o(1)},\qquad
 |G_v/G_d|\le C_A\sqrt{L/N}.                                        \tag{11}
\]

Indeed the base response may first be restricted to `C_D` at a loss smaller
than any power of `N`; (6) then compares it from below with the shifted
response.  The central estimates in (1) similarly give, with
`beta=G_{ad}/G_d=O_A(L)`,

\[
 |G_{av}-\beta G_v|\le C_A L^{3/2}G_d/\sqrt N.                     \tag{12}
\]

The target Ward/residual decomposition is exact:

\[
 G_{aa}={N\over2}m''(a)G_d+E_{a,u}[H_0 Z_a].                        \tag{13}
\]

Here `E_{a,u}` is expectation under the shifted target product law and
`E_{a,u}Z_a=0`.  Consequently the last expectation is exactly
`E_{a,u}[(H_0-1)Z_a]`.  On `C_D`, (6) transfers the source `L^2` profile
estimate to the target law, while the r054 target-law estimate gives
`||Z_a||_2=O(L^2)`.  On `C_D^c`, use the ordinary target Chernoff and
Gamma-moment estimate (9), not the response-tail lemma.  Hence

\[
 |E_{a,u}[(H_0-1)Z_a]|\le N^{-1/4+o(1)}L^2=o(NG_d).
\]

Since `m''(a)=-p_mu+o(1)`, (13) yields

\[
 -C_A NG_d\le G_{aa}\le-c_A N G_d.                                  \tag{14}
\]

Implicit differentiation now gives

\[
 d'=-G_v/G_d=O_A(\sqrt{L/N}),\qquad
 a'=-{G_{av}-\beta G_v\over G_{aa}}
       =O_A(L^{3/2}/N^{3/2}).                                       \tag{15}
\]

Choose `B>2A\max\{C_A,1\}`.  Integrating (15) gives strictly smaller bounds
than (10) for large `N`, so the maximal-open-interval argument continues the
unique even regular solution over the entire natural window
`|v|<=A sqrt N`.

Finally let `d(a,v)` be the solution of `G(d,v,a)=1`, and put `y=d_a`.
The target location derivative gives the exact response identity

\[
 G_{dd}=-B_R(k_{\rm source,+},\widetilde k_+-2),\qquad
 \gamma:={G_{dd}\over G_d}=O_A(\sqrt{L/N}),                         \tag{16}
\]

where (16) uses exactly the last weight in (5).  With
`\beta=G_{ad}/G_d=O_A(L)`, a second implicit differentiation gives the
Riccati identity

\[
 y'=d_{aa}=-{G_{aa}\over G_d}-2\beta y-\gamma y^2.                 \tag{17}
\]

At the stationary receiver `a_*=a(v)`, `y(a_*)=0`.  On the interval
`|a-a_*|<=C_A L^{3/2}/N`, bootstrap (17) first gives
`|y|<=C_A L^{3/2}`.  The two lower-order terms in (17) are then
`O_A(L^{5/2})+O_A(L^{7/2}/sqrt N)=o(N)`.  By (14), for large `N`,

\[
 0<c_A N\le d_{aa}\le C_A N.                                       \tag{18}
\]

Thus `a_*` is the local minimizer of `d(a,v)`.  Taylor's theorem and (15)
give

\[
 0\le d_{\rm fix}(v)-d_{\rm opt}(v)
 \le C_A N|a(v)-c|^2=O_A(L^3/N)=o(1).                                \tag{19}
\]

Since the Pareto-cone volume is an integral plus endpoint values of
`e^{-2d}`, the uniform `o(1)` difference in (19) yields a volume ratio
`1+o(1)` on this central arc.  This is a stability statement only; it does
not supply a nonzero lower bound for the gain.
