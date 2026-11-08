# A broad radial window closes the bounded-input screening module

## Statement and scope

Fix a cap \(0\le F\le A\), a reference action
\[
 c=\mu+O(\sqrt{\log N/N}),
\]
and a broad interval
\[
 W_N=[N\mu-C\sqrt{N\log N},\,N\mu+C\sqrt{N\log N}],
\]
where \(C\) may be chosen after a prescribed tail power.  Define the
physical local cost
\[
 j_WG(u,X,\sigma)
 =\int {\bf1}_{W_N}(s)e^{-s}
 [K_c^sG_{u+s}](X,\sigma)\,ds,\qquad G=\Phi_N(F).
\]
This note verifies a bounded-input *screening and scalar-window-capture*
module.  It does not provide the missing common-source comparison,
reference-to-moving-action step, anchor geometry, or halo packing needed
for a general upper bound.

## Uniform radial tail and mean repair

For \(c\) in the displayed compact shrinking window, the increment
exponential moments are uniformly bounded in a fixed neighborhood of zero.
Its mean differs from \(N\mu\) by \(O(\log N)\), and its variance is
\(\asymp N\).  Uniform Bernstein--Chernoff bounds therefore give, for
each prescribed \(M\), a choice of \(C=C(M)\) such that
\[
 \sup_c\mathbb P_c\{S\notin W_N\}\le N^{-M}.                \tag{1}
\]
This applies to the full Gamma--Bernoulli radial density
\(p_c(s)=e^{-s}H_c(s)/Z_c\), not merely to a Gaussian proxy.

At a receiver, put
\[
 h(s)=[K_c^sF_{u+s}](X,\sigma),\qquad h_W=h{\bf1}_{W_N}.
\]
If the actual receiver average is at least one, (1) and \(F\le A\) give
\[
 \int p_c(s)h_W(s)\,ds\ge1-\delta_N,\qquad
 \delta_N=A N^{-M}.                                        \tag{2}
\]
Take \(M>1\) and set \(\widetilde h=h_W/(1-\delta_N)\).  Then
\(\int p_c\widetilde h\ge1\), \(\widetilde h\le A+1\) for large \(N\),
and the elementary scaling inequality
\[
 \Phi_N(ax)\le a^N\Phi_N(x)\qquad(a\ge1)
\]
gives
\[
 \int e^{-s}\Phi_N(\widetilde h(s))\,ds
\le(1-\delta_N)^{-N}\int_{W_N}e^{-s}\Phi_N(h(s))\,ds
=(1+o(1))\int_{W_N}e^{-s}\Phi_N(h(s))\,ds.                 \tag{3}
\]

The last cost is controlled by the original common source, not by a
receiver-dependent modification of it.  Positivity and Jensen give
\[
 \Phi_N(h(s))\le[K_c^s\Phi_N(F_{u+s})](X,\sigma)
 =[K_c^sG_{u+s}](X,\sigma),
\]
hence
\[
 \int_{W_N}e^{-s}\Phi_N(h(s))\,ds\le j_WG(u,X,\sigma).      \tag{4}
\]
Thus the receiver-dependent \(\widetilde h\) is only an auxiliary scalar
profile; every cost bound in (3) is paid by the same unmodified \(G\).

## Physical screening cost

For physical measure \(d\mathfrak m=e^{-u}du\,dH\,d\lambda\), the
individual \(s\)-term in \(j_W\) is \(L^1\)-isometric after the
\(u\mapsto u+s\) change of variables.  Consequently
\[
 \|j_WG\|_{L^1(d\mathfrak m)}
\le |W_N|\,J(G)
=N^{1/2+o(1)}J(G).                                         \tag{5}
\]
At threshold \(t=J_{\rm ref}N^\varepsilon\), Markov gives
\[
 \mathfrak m\{j_WG>t\}
\le N^{1/2-\varepsilon+o(1)}\,{J(G)\over J_{\rm ref}}.     \tag{6}
\]
This is a genuine physical weak screening estimate with exponent
\(\gamma=1/2\), obtained only because the radial integration has been
restricted to \(W_N\).

## Compatibility with r069 capture

The critical forbidden band \(I\) of r069 lies inside \(W_N\) for all
large \(N\).  Apply its cap-\((A+1)\) capture statement to
\(\widetilde h\).  Equations (3)--(4) verify its global scalar-budget
hypothesis whenever \(j_WG\le J_{\rm ref}N^\varepsilon\); no uncontrolled
tail cost remains because \(\widetilde h\) was set to zero outside \(W_N\).
The harmless factor in (3) may be absorbed into \(N^{o(1)}\).

Since
\[
 {1\over2}<2\alpha={9-\mu\over12},
\]
one may choose
\[
 K>\max\{\alpha,1/4\},\qquad {1\over2}<\varepsilon<2K.
\]
For the near-optimal shifted reference take \(K=\alpha+\eta\) with
\(\eta>0\) arbitrarily small.  The screened bad set in (6) then decays,
and the r069 retained-band loss is only \(N^{\eta+o(1)}\).  This closes
the bounded-input local-cost screening/window-capture module without a
full-\(j\) weak-\((1,1)\) theorem.

What remains outside this module is substantial: translating the captured
scalar profile into a common physical source comparison, controlling
reference versus moving actions, and packing the resulting receiver
anchors.  Therefore this result must not be stated as a new general
\(C_N\) upper bound.
