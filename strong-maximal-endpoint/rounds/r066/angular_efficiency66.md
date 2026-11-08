# Verified angular efficiency lemma and the remaining global step

## Outcome

The usable result is a fixed-radial, finite-depth, labelled gain--Bregman
inequality. It charges a policy's conditional angular gain to an angular
Bregman defect, retaining both the exact Dirichlet Jacobian and the
likelihood envelope from `r064/envelope64.md`.

It does not turn near-unit receiving means plus a global critical modular
budget into the action-independent terminal required in
`r065/universal_action65.md`. The exact missing implication is Section 3.

Let `d_j=b_j-1`, fix receiver radii `V_j>0` and block sums `s`, and write
`w=u+s` for the source radial coordinate. Orthant labels are retained and
are averaged with their uniform probability.

## 1. Fixed-s gain--Bregman inequality

Let `K_c^{V,s}` be the labelled conditional angular Markov kernel. Its
transport from receiver to source angles is

\[
 \eta_j=\frac{V_j\xi_j+z_j}{V_j+s_j}.                   \tag{1}
\]

For an action-independent radial baseline `B=B(w)`, choose an affine
supporting function, not merely a scalar subgradient,

\[
 \ell_B(x)=\Phi_N(B)+g_B(x-B),\qquad g_B\in\partial\Phi_N(B),
 \quad D_B(x)=\Phi_N(x)-\ell_B(x)\ge0.                  \tag{2}
\]

Put `d_B(epsilon)=D_B(B+epsilon)`. Since `ell_B` is affine and the kernel
is Markov, conditional Jensen gives

\[
 D_B(K_c^{V,s}F(\xi))\le K_c^{V,s}(D_B(F))(\xi).         \tag{3}
\]

For a measurable action policy `pi(xi,sigma)` in the near-`mu` box,

\[
\begin{aligned}
 &d_B(\varepsilon)\mathbb E_\sigma\,\sigma_\angle
 \{\xi:K_{\pi(\xi,\sigma)}^{V,s}F_\sigma(\xi)
                  \ge B+\varepsilon\}\\
 &\quad\le \mathbb E_\sigma\int K_{\pi(\xi,\sigma)}^{V,s}
 D_B(F_\sigma)(\xi)\,d\sigma_\angle(\xi)\\
 &\quad\le
 \left\{\prod_j\left(\frac{V_j+s_j}{V_j}\right)^{d_j}\right\}
 \widetilde Z_s\mathbb E_\tau\int D_B(F_\tau)(\eta)
 \,d\sigma_\angle(\eta).                               \tag{4}
\end{aligned}
\]

Here `sigma_angle` is uniform probability on the product simplices and

\[
 \widetilde Z_s\le\prod_j[1+b_j\log(\beta_j/\alpha_j)].           \tag{5}
\]

Here [\(\alpha_j,\beta_j\)] is the action interval, distinct from the angular dimension d_j. The last factor is the joint-sign envelope in (27)--(28) of
`envelope64.md`. To prove (4), first use (3), condition on increments and
signs, and remove the policy with that envelope. For each fixed increment,
the exact Dirichlet push has the displayed factor
`prod_j((V_j+s_j)/V_j)^{d_j}`; enlarge the shifted source simplex to the
full source simplex. Finally the map `sigma -> sigma(-1)^J` permutes
labels. This proves (4), including the finite-depth factor.

Multiplication by the receiver radial Jacobian `prod_j V_j^{d_j}` cancels
this factor against the source Jacobian `prod_j(V_j+s_j)^{d_j}`, exactly as
in `angular_budget64.md`. That cancellation requires radial weighting; it
is not a literal unweighted angular contraction at fixed `V,s`.

At `B=1`, take the left supporting slope `g_B=1`. For `epsilon=t/N`,

\[
 d_1(t/N)=\Phi_N(1+t/N)-1-t/N=e^t-1+o(1)                 \tag{6}
\]

for fixed `t`. This is the exact exponential gain--cost scale.

## 2. Local near-efficient action dependence

For an angular set `A` of probability `m`, let

\[
 b_+=t/N,\qquad a_-=\frac{mb_+}{1-m},\qquad
 R_A=(1-a_-)1_{A^c}+(1+b_+)1_A.                           \tag{7}
\]

It has angular mean one and

\[
 \int D_1(R_A)d\sigma
 =m\{\Phi_N(1+t/N)-1-t/N\}=m\{e^t-1+o(1)\}.             \tag{8}
\]

The conditional terminal is

\[
 K_c^{V,s}R_A(\xi)=1-a_-+(a_-+b_+)K_c^{V,s}1_A(\xi).     \tag{9}
\]

It is action-dependent whenever two conditional kernels differ on `A`.
This already happens in one block because

\[
 q_c(z)=\frac{e^{-z}(1+1_{z>c})}{1+e^{-c}}               \tag{10}
\]

varies with `c`; a deep interior chart supplies an open separating set.
Thus near-unitness cannot establish an exact action-independent terminal
identity.

For fixed `t` and `m -> 0`, the deep physical product

\[
 F(w,\eta)=\mathbf1_{w\in\mathbb R_+^r}h_0(\sum_jw_j)R_A(\eta)                       \tag{11}
\]

has `J(F)=(1+O(m))J(h_0)` after coarea (use any fixed cap above two, or
truncate the negligible top of the ramp), and its receiving means differ
from those of `h_0` by `O(1/N)`. This is only a local obstruction to an
identity. It does not refute an aggregate common-action reduction. In
particular `m=1/\widetilde Z_s` gives defect `O(1/\widetilde Z_s)`, which
is compatible with a condition of size `N^{-r/2+o(1)}`.

## 3. Exact missing implication

The r065 theorem concerns receiver averages over the random sums, whereas
(4) is fixed-`s`. Integrating its conditional exceptional sets requires
the policy-dependent sum law; moreover a small exceptional angular set at
each `s` need not give a negligible physical halo after receiver-cone
reuse. Neither conclusion follows from (4) by Fubini alone.

What is missing is a joint assertion with three compatible parts: a single action-independent block profile whose translates describe all retained receivers, a comparison of its actual global modular with the physical source modular, and a cone allocation that makes the discarded physical receiver volume negligible relative to Q_N J(F). A fixed-s exceptional-angle bound does not supply these parts. In particular, separately choosing a baseline B_u for each receiver would not give the one global source profile needed by r065.

The spatial repair of r065 can absorb an o(log N/N) received-mean loss only after those representation and budget conditions have been proved. The present gain--Bregman inequality is an interface for that task, not a proof of the missing representation or a general C_N upper bound.
