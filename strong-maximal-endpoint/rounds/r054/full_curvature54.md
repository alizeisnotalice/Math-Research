# A scoped full-quadrant contrast-curvature theorem

This proof concerns only the capped, two-equal-block, one-receiver convex
budget below.  It does **not** imply a bound for the strong maximal constant
`C_N` or for a multi-receiver/variable-rectangle problem.

## Theorem

Let `N=2b`, `b -> infinity`, and let
\[
 p_c=(1+e^c)^{-1},\quad m(c)=1+c p_c,
 \quad \mu=1+e^{-\mu},\quad p_\mu=(\mu-1)/\mu.
\]
Let `f_{b,c}` be the density of the sum of `b` independent copies of
`E+cB_c`, where `E` is unit exponential and `B_c` is Bernoulli(`p_c`), and
put
\[
 R_{b,c}(t)=e^t f_{b,c}(t),\qquad
 d\xi=e^{-t_1-t_2}dt_1dt_2,
 \qquad R(t_1,t_2)=R_{b,c}(t_1)R_{b,c}(t_2).
\]
For
\[
 \Phi_N(h)=h(1+\log^+h)^{N-1},\qquad 0\le h\le2,
\]
let `P(c_1,c_2)` be the actual full-positive-quadrant cap-two budget at
target one,
\[
 \inf\left\{\int\Phi_N(H)d\xi:0\le H\le2,
 \ \int H R_{b,c_1}(t_1)R_{b,c_2}(t_2)d\xi\ge1\right\}.
\]
Assume the standard local KKT regularity at the symmetric point
`(c,c)` (positive locally unique multiplier, regular lower level, and the
permitted differentiation/coarea interchanges).  Uniformly when
\[
 |c-\mu|\le C\sqrt{\log b/b},\tag{1}
\]
the contrast path `theta -> (c+theta,c-theta)` satisfies
\[
 \boxed{\quad (\log P)''_-(c,c)=N p_\mu+o(N).\quad}\tag{2}

The fixed finite exceptional values of `b` are irrelevant to this
asymptotic statement.

## 1. Exact KKT--Ward reduction

Let `lambda>0` be the KKT multiplier, `a=lambda R`, and `H=h(a)` the
cap-two response:
`H=0` for `a<1`, `H=1` for `1<a<=N`, the inverse of `Phi_N'` in the soft
region, and `H=2` at the cap.  Write `nu=R dxi`,
`s=partial_c log R_{b,c}`, and `k=partial_t log R_{b,c}`.  The Fenchel
density `Phi_N(H)-lambda RH` has no delta mass at the lower jump and has
vanishing axis/infinity boundary contributions on the positive quadrant.
Integration by parts gives the Ward identity
\[
 \lambda E_\nu H(1-k)=P.\tag{3}
\]

For one block define
\[
 e(t)={R_{cc}(t)\over R(t)}-bm''(c)(1-k(t)),
\]
and for two independent blocks set
\[
 Z=e(T_1)+e(T_2)-2s(T_1)s(T_2).
\]
The exact symmetric contrast Hessian calculation, including both the lower
coarea response and the soft inverse-curvature response, is
\[
 (\log P)''_-=N[-m''(c)]-{\lambda\over P}(\mathcal E_N+B_{--}),
 \qquad \mathcal E_N=E_\nu[HZ].\tag{4}
\]
Thus it suffices to prove
\[
 {\lambda\over P}B_{--}=o(N),\qquad
 {\lambda\over P}|\mathcal E_N|=o(N).\tag{5}

## 2. Compact Fourier bridge and a drifting Fisher estimate

Let `g_{n,d}` be the density of the centered sum of `n` copies of
`E+dB_d`.  On every compact `D` in `(0,infinity)`, uniformly for bounded
`x`, Fourier inversion gives the differentiated local limit estimate
\[
 g_{n,d}(x)=(2\pi n v(d))^{-1/2}\{1+O_{C^1(D,x)}(n^{-1})\}.
\tag{6}
\]
Here `C^1` includes one `d` derivative and one bounded-`x` derivative.
For completeness, let `psi_d` be the centered characteristic function and
put `zeta=u/sqrt n`.  On `|u|<=n^(1/12)`,
\[
 e^{-i\zeta x}\psi_d(\zeta)^n=e^{-v(d)u^2/2}
 \left[1-{i(xu+\kappa_3(d)u^3/6)\over\sqrt n}
 +O_D((1+|u|^6)/n)\right].\tag{7}
\]
The first correction is odd imaginary.  Parameter differentiation changes
the leading Gaussian amplitude by `-v'(d)u^2/2`; after that leading term is
accounted for, the same odd cancellation proves the relative `O(1/n)`
remainder.  An `x` derivative supplies `-iv`, whose leading integral is
also odd.  The middle frequencies are exponentially small from
`|psi_d(v)|<=exp(-c_Dv^2)`.  At outer frequencies,
`|psi_d(v)|<=(1+v^2)^(-1/2)`; compact positivity of the Bernoulli factor
gives `|partial_d log psi_d(v)|<=C_D(1+|v|)`, so the differentiated outer
tail is still exponentially small after four powers are retained for
integrability.  This proves (6).

The exact common-tilt score representation is
\[
 s_{b,c}(t)=bp_c\left[1-{u(d)B_b(d)\over u(c)}\right],
 \quad {t\over b}={c m(d)\over d},\quad
 B_b(d)={g_{b-1,d}(m(d)-d)\over g_{b,d}(0)}.\tag{8}
\]
By (6), `B_b=1+O_{C^1}(b^{-1})` on compact `d` sets.  Since `u` has a
strict quadratic maximum at `mu`, the local inverse in (8), the global
lower score estimate, and the elementary upper bound `s<=bp_c` imply
\[
 |s_{b,c}(t)|\le C\left(1+b(c-\mu)^2+
 {\{t-bm(c)\}^2\over b}\right).\tag{9}
\]
The estimate first holds on a fixed compact `t/b` interval; outside it the
last term is `>=c b`, so the global one-sided score bound supplies the same
conclusion.  Uniform fourth moments of the increments then yield the
drifting Fisher bound
\[
 I_b(c):=E s_{b,c}(T)^2\le C(1+b(c-\mu)^2)^2=O((\log N)^2).\tag{10}

## 3. A global multiplier estimate

For the cap-two KKT solution,
\[
 {\lambda\over P}\le2N\tag{11}
\]
for all large `N`, uniformly in (1).  Here is a short proof.  Off the cap,
`aH<=N Phi_N(H)`, whence
\[
 \lambda\{1-2\nu(H=2)\}\le NP.\tag{12}
\]
Take the fixed compact interval `J=[6/5,27/20]`, which contains (1) for
large `b`.  The `k=0` Gamma summand gives uniformly on `J`
\[
 \log R_{b,c}(6b/5)>.90b.\tag{13}
\]
Both blocks exceed `6b/5` with probability tending to one exponentially
fast.  If `lambda e^(.90N)>=Phi_N'(11/10)`, then `H>=11/10` on an event of
probability above `.95`, contradicting `E_nu H=1`.  Hence that inequality
is false.  On the cap,
\[
 \log R\ge .90N+\log{\Phi_N'(2)\over\Phi_N'(11/10)}>1.32N.\tag{14}
\]
Each increment density is at most one, so `f_{b,c}<=1` and (14) implies
`T_1+T_2>1.32N`.  Since `m(c)<=mu<1.3`, a uniform Chernoff bound makes the
cap probability exponentially small.  Equation (12) gives (11).

## 4. The entire response term is negligible

On the soft set,
\[
 {\Phi_N'(h)\over\Phi_N''(h)}
 ={h(1+\log h)(N+\log h)\over
 (N-1)(N-1+\log h)}\le{6\over N},\qquad1\le h\le2.\tag{15}
\]
Block independence and `Es=0` therefore give
\[
 B^{\rm soft}_{--}\le{12\over N}I_b(c)=O((\log N)^2/N).\tag{16}

For the lower jump, set `U=log R_{b,c}(T)` and let `q_b` be its density.
The positive truncated-polynomial representation gives
`k_b(t)>=(b-1)/t`.  Since the block density is a binomial mixture of
shifted Gamma densities,
\[
 \sup_t t f_{b,c}(t)=O(\sqrt b),\qquad \|q_b\|_\infty\le C/\sqrt b.\tag{17}
\]
Indeed `x gamma_b(x)=b gamma_{b+1}(x)` and
`sup gamma_b=O(b^{-1/2})`, the latter following directly from the Gamma
characteristic-function integral.  At the lower jump,
`delta(lambda e^{u_1+u_2}-1)=delta(log lambda+u_1+u_2)`: there is no extra
Jacobian.  Thus a convolution bound gives, uniformly in the level,
\[
 B^{\rm low}_{--}\le4\|q_b\|_\infty I_b(c)
 =O((\log N)^2/\sqrt N).\tag{18}

Equations (11), (16), and (18) prove the first part of (5).

## 5. The KKT profile is close to one in `L^2(nu)`

The density of one increment is within a factor two of the unit exponential
density.  Comparison with the exponential Poincare inequality gives
`Var F(Y_c)<=8 E|F'(Y_c)|^2`.  Tensorization and
\[
 k(T)-1=E\left[{b-1\over G_b}-1\mid T\right],\qquad
 E\left({b-1\over G_b}-1\right)^2={1\over b-2},\tag{19}
\]
give, for `W=log R_b(T_1)+log R_b(T_2)`,
\[
 \operatorname{Var}W\le12N.\tag{20}

Let `q=-log lambda` and `g(u)=(h(e^u)-1)_+`.  Directly from (15),
\[
 \min\{(u-\log N)_+/N,1\}\le g(u)
 \le6(u-\log N)_+/N.\tag{21}
\]
If `M=EW` and `a=N^(3/4)`, the identities
\[
 P_0:=P(W<q)=Eg(W-q),\qquad E|H-1|=2P_0\tag{22}
\]
follow from `EH=1`.  If `q<=M-a`, Chebyshev makes `P_0=O(N^{-1/2})`, whereas
on `W>=M-a/2`, (21) gives `g(W-q)>=(a/2-log N)/N>=a/(3N)` eventually, a
contradiction.  Hence `q>=M-a`; using (20)--(22),
\[
 P_0\le {6\over N}E(|W-M|+a)=O(N^{-1/4}),
 \qquad \|H-1\|_2=O(N^{-1/8}).\tag{23}

## 6. A logarithmic `L^2` bound for the kernel residual

At fixed `t`, implicit differentiation in (8) gives
\[
 d_c={dm(d)\over cv(d)},\qquad
 tk+cs=b m(c)-1.\tag{24}
\]
Writing `F(c,d)=p_c-u(d)/c`, define
\[
 G=F_c+F_d{dm(d)\over cv(d)}-m''(c)(1-d/c).
\]
Direct algebra gives
\[
 G(c,c)=-{m'(c)^2\over v(c)},\qquad
 G(\mu,\mu)=G_c(\mu,\mu)=G_d(\mu,\mu)=0.\tag{25}
\]
Also `F(mu,mu)=F_c(mu,mu)=F_d(mu,mu)=0`.  The bridge estimate in §2 therefore
implies, on
`|T-bm(c)|<=sqrt(Cb log b)`,
\[
 s=O(\log b),\qquad e=s^2+bG+O(1)=O((\log b)^2).\tag{26}

For the complement, use the exact fixed-`(t,K)` identity under
`T=G_b+cK`,
\[
 {R_{cc}(T)\over R(T)}=E\left[
 b^2p_c^2-bp_c(1-p_c)-2bp_c(b-1){K\over G_b}
 +(b-1)(b-2){K^2\over G_b^2}\ \middle|\ T\right].\tag{27}
\]
The inverse Gamma moment `E G_b^{-8}=O(b^{-8})` gives
`||R_{cc}/R||_4=O(b^2)`.  The same mixture formula gives `||k||_4=O(1)`;
together with the global score bounds this proves `||Z||_4=O(b^2)`.
Uniform exponential moments give an arbitrarily high polynomial tail for
the complement of the central band.  Hölder plus (26) therefore yields
\[
 \|Z\|_2=O((\log b)^2).\tag{28}

Finally `EZ=0`: normalization gives `E(R_{cc}/R)=0`, equation (19) gives
`Ek=1`, and the cross-score term factors.  Thus (23) and (28) imply
\[
 |\mathcal E_N|=|E(H-1)Z|
 =O(N^{-1/8}(\log N)^2).\tag{29}

By (11), this proves the second part of (5).  Since
`-m''(c)=p_mu+o(1)` under (1), substitution into (4) proves (2).

## Scope and verification

The proof uses the actual full positive quadrant, cap `A=2`, a single
receiver, two equal independent blocks, and the stated regular KKT
hypotheses.  It is not a theorem about `C_N`.

The algebraic central cancellations in (25) were independently checked by
the actual Wolfram 15 script
`code/z_central_exact54.wl`, which returns
the diagonal identity and the two first derivatives as exact zero after
the defining relation for `mu` is substituted.  This computation is
corroboration only; the displayed proof is analytic.
