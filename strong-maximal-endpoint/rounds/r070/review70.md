# Review: the weak-screen exponent required by band capture

## Conclusion

The proposed relaxation is algebraically correct.  One does **not** need a
dimension-free weak `(1,1)` estimate for the local cost operator.  At a fixed
reference action, a weak bound with constant

\[
 \|\mathcal J_c\|_{L^1\to L^{1,\infty}}
 \le W_N=N^{\gamma+o(1)},\qquad \gamma\le2\kappa,
 \tag{W}
\]

is enough at the leading exponent for the efficiency-qualified screening plus
the r069 band-capture argument.  The strict version `\gamma<2\kappa` is a
convenient sufficient form with a fixed gap.  This is a conditional
implication: the present elementary `L^p` argument gives only an `O(N)`
*threshold-specific capped screening coefficient*, and neither proves (W)
nor proves a matching lower bound for an actual uniform weak constant.

The simple local numerical check gives

\[
 \alpha=(9-\mu)/24=0.32173064405162194,\qquad
 2\alpha=0.6434612881032439.
\]

Thus a genuinely proved `O(\sqrt N\log N)` weak constant would have
`\gamma=1/2<2\alpha`; it would suffice for the shifted fixed reference.
For the unshifted reference the corresponding threshold is
`2\beta=(3-\mu)/3\approx0.5738451524`, so the same hypothetical bound
would also suffice there.

## Screening calculation

Write `\nu` for the physical measure (including the normalized angular-label
Haar factor), `G=\Phi_N(F)`, and `\mathscr J=\int G\,d\nu`.  Suppose the
terminal cap gives `0\le F\le2`, so that

\[
 0\le G\le M_N:=\Phi_N(2)=2(1+\log2)^{N-1}.
\]

At the chosen *fixed* action, let

\[
 j=\mathcal J_cG,\qquad t=J_{*,2}N^\varepsilon.
\]

From (W),

\[
 \nu\{j>t\}
 \le {W_N\mathscr J\over t}
 =Q_{\rm ref}N^{\gamma-\varepsilon+o(1)}\mathscr J.
 \tag{1}
\]

On the complement, `j\le J_{*,2}N^\varepsilon`, precisely the extra local
budget hypothesis needed in the r069 scalar capture theorem.  Choose

\[
 K=\kappa+\delta,\qquad \gamma<\varepsilon<2K.
 \tag{2}
\]

Then r069 fixed-action transport controls the successful good set by

\[
 \nu(E_{\rm good})
 \le Q_{\rm ref}N^{K-\kappa+o(1)}\mathscr J
 =Q_{\rm ref}N^{\delta+o(1)}\mathscr J,
 \tag{3}
\]

while (1) gives a strictly smaller bad-set contribution.  If
`\gamma<2\kappa`, one can first choose a sufficiently small fixed
`\delta>0` and then choose `\varepsilon` as in (2).  At the endpoint
`\gamma=2\kappa`, take `K=\kappa+\delta` and
`\varepsilon=2\kappa+\delta`.  Then the bad-set exponent is `-\delta`,
while the good-set loss is `+\delta`.  Taking `N\to\infty` first and then
`\delta\downarrow0` yields the same leading `\kappa` exponent.  Thus
`\gamma\le2\kappa` is the leading threshold; strict inequality merely
avoids this limiting bookkeeping.

## What the available `L^p` estimate actually says

The known Markov/Minkowski estimate is

\[
 \|\mathcal J_cG\|_p\le p'\|G\|_p\qquad(p>1).
\]

It and the cap only imply

\[
 \nu\{\mathcal J_cG>t\}
 \le (p')^p\,(M_N/t)^{p-1}\,{\mathscr J\over t}. \tag{4}
\]

For `H=\log(M_N/t)\asymp N`, the convenient choice
`p=1+1/H` makes the coefficient in (4)

\[
 e(H+1)^{1+1/H}=O(N).
\]

Optimizing in `p` does not change its order: if `q=p-1`, the logarithm of
the coefficient is

\[
 (1+q)\log((1+q)/q)+qH,
\]

whose minimizer has `q\asymp1/H` and minimum of order `\log H`.
Consequently this method supplies a capped, threshold-specific *upper
screening coefficient* of order `N`, not a lower bound (or a uniform upper
bound) for the actual weak constant.  It therefore misses both shifted and
unshifted thresholds (`2\kappa<1`) but creates no obstruction to a better
geometric proof.

## Scope that remains necessary

1. The weak estimate must apply to the actual fixed finite reference action
   and the full label space.  A result only for the `c=\infty` uniform-simplex
   model does not by itself screen the shifted finite-`c_N` reference.
2. The result only screens receivers for which the fixed-reference formulation
   is already legitimate.  It does not provide the still-missing
   variable-action-to-common-reference transfer, its `o(\log N/N)` mean
   accuracy, or a compatible moving-anchor allocation.
3. This uses the weak estimate for `j=\mathcal J_c\Phi_N(F)`, not a radial
   band maximal envelope and not an assertion that global efficiency alone
   gives a pointwise local budget.  The rogue-packet qualification counterexample
   in r069 is unaffected.
