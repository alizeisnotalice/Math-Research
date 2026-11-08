# r075 review: coarse global screen for arbitrary critical vector policies

## Verdict

The coarse lemma is valid, provided the wide band and its shifted copy stay
inside a fixed positive macroscopic sum compact set.  It is precisely a
high-threshold screen: it charges vector policies with exponentially large
windowed cost without any differentiability, injectivity, or common-anchor
assumption.  Its loss is `exp(O(sqrt(N log N)))`, so it cannot replace the
polynomial-slack (`N^epsilon`) screen needed for the sharp core.

## Conditional-kernel domination

Put `c_-=c_r-Krho`, `c_+=c_r+Krho`, with
`rho=sqrt(log(N)/N)`, and let every coordinate of a measurable vector action
`mathbf c` lie in `[c_-,c_+]`.  For an increment-label set `J`, use the
unnormalized conditional simplex domain

\[
 D_{\mathbf c,J}(s)=\{z_i\ge c_i\mathbf1_{i\in J}:\ \sum_i z_i=s\}.
\]

The conditional-label cancellation makes its unnormalized density just
simplex measure: the Bernoulli label factors cancel the exponential shift.
Consequently, label by label,

\[
 D_{c_+,J}(s)\subset D_{\mathbf c,J}(s)\subset D_{c_-,J}(s),
 \qquad H_{\mathbf c}(s)\ge H_{c_+}(s).                  \tag{1}
\]

The second inclusion bounds the unnormalized numerator of a positive
conditional convolution by the `c_-` numerator.  Dividing by the first
lower bound in (1) gives the correctly directed operator inequality

\[
 0\le K_{\mathbf c}^{,s}G
 \le {H_{c_-}(s)\over H_{c_+}(s)}K_{c_-}^{,s}G.           \tag{2}
\]

This also holds for the complete joint label convolution: the comparison is
made separately for each increment label and label translation is a
permutation of the finite output label group.  Thus there is no missing
`2^N` factor.

Let `Delta=c_+-c_-=2Krho`.  Termwise in the simplex sum,

\[
 H_{c_+}(s)\ge H_{c_-}(s-N\Delta).                        \tag{3}
\]

Indeed, for a term with `k` jumps, the difference between its two positive
parts is `(N-k)Delta>=0`.  On the expanded macroscopic compact containing
`W_N` and `W_N-NDelta`, the established radial logarithmic-rate bound gives
`partial_s log H_{c_-}(s)<=B`.  Hence

\[
 E_N:=\sup_{s\in W_N}{H_{c_-}(s)\over H_{c_+}(s)}
 \le\exp\{BN\Delta\}=
 \exp\{O(\sqrt{N\log N})\}.                               \tag{4}
\]

No smoothness of the policy enters (1)--(4).

## Physical Fubini estimate

For the policy windowed cost, (2) gives pointwise

\[
 j_{W_N}^{\rm policy}G(u,X,\sigma)
 \le E_N\int_{W_N}e^{-s}
 [K_{c_-}^sG_{u+s}](X,\sigma)\,ds.                        \tag{5}
\]

Integrate against the physical measure `e^{-u}du dX dlambda`.  With the
source coordinate `w=u+s`, the two exponential weights combine exactly:

\[
 e^{-u}e^{-s}=e^{-w}.
\]

The fixed `c_-` complete joint kernel is Markov on `X,sigma`, so Tonelli
gives

\[
 \int j_{W_N}^{\rm policy}G\,d\nu
 \le |W_N|E_N\int G\,d\nu.                               \tag{6}
\]

This is the required common-global-anchor bound; it does not use a receiver
change of variables.  In particular it is immune to the folding and
pullback failures of re-anchoring.

At a threshold `J_ref exp(aN)`, Markov and (6) give the advertised form

\[
 \nu\{j_{W_N}^{\rm policy}G>J_{\rm ref}e^{aN}\}
 \le Q_{\rm ref}\,e^{-aN+o(N)}\,\mathcal J(G),           \tag{7}
\]

where the `o(N)` absorbs `log|W_N|+log E_N=O(sqrt(N log N))`.
Therefore a construction such as the optimizer example, whose windowed
cost is `exp(nu N+O(sqrt(N log N)))J_ref`, is globally chargeable by choosing
any fixed `0<a<nu`.

## Limitation

For a threshold only `J_ref N^epsilon`, (6) yields the factor
`exp(O(sqrt(N log N)))N^{-epsilon}`, which is not small.  The coarse lemma
therefore does not supply the polynomial-slack local-cost qualification
needed in the sharp argument.  It corrects the earlier overstatement that
the r072 fixed-action screen alone removed arbitrary variable-action
optimizers, while leaving the general sharp-core/action-packing problem
open.
