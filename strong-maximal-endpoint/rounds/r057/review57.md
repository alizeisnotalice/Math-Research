# r057 independent audit — common-window response normalization

Date: 2026-11-07.  Scope is the regular, capped, equal two-block,
single-receiver KKT model only.

Let `r=lambda/P`, `g=k_1+k_2-2`, and let `B` be the full lower-level plus
soft response form.  The Ward calculation gives exactly
\[
 B(1,k_1+k_2)=2/r,
 \qquad B_{00}=1/r-\tfrac12B(1,g).
\]
For the lower term, disintegrate with respect to the two block log-kernel
values.  Conditional Cauchy (so no injectivity of `log R_b` is being
assumed) gives
\[
 B_{\rm low}(g,g)\le4\|q_b\|_\infty E(g^2).
\]
The stated `q_b=O(N^{-1/2})` and location Fisher `E(g^2)=O(N^{-1})` make
this `O(N^{-3/2})`; the soft term is `O(N^{-2})`.  Hence
\[
 |B_{00}-1/r|\le C\sqrt{B_{00}}N^{-3/4}.
\]
Since `r<=2N`, setting `y=sqrt(r B00)` implies
`|y^2-1|<=C y N^{-1/4}`, and therefore
\[
 rB_{00}=1+O(N^{-1/4}).
\]
This normalization algebra is valid and is the key new ingredient.

With `M=O(N^{-1/2}(log N)^{3/2})` and
`B_{++}=O((log N)^2/sqrt N)`, the three nonstationary Schur terms are,
respectively,
\[
 O(N^{3/4}(\log N)^3),\quad
 O(N^{3/4}(\log N)^{5/2}),\quad
 O(\sqrt N(\log N)^2),
\]
and hence are all `o(N)`.  Wolfram 15.0 was actually run to simplify these
three powers; it returned exactly these exponents.  Subject to the existing
uniform `Z_+` and response bounds, this closes the off-stationary obstruction.

One presentation requirement for the profile-mass lemma: keep distinct
\(\pi_0=\Pr(W<\ell)\) and \(E g(W-\ell)\).  Their equality is a KKT
mean identity, not a definition.  The correct contradiction chain is:
an overly low `ell` first gives `Pr(W<ell)<=N^-M` from the moderate lower
bound for `W`; the mean identity then gives `E g<=N^-M`; the central lower
bound on `W-ell` makes `g` of order `sqrt(log N/N)` on a fixed-positive
probability event, contradiction.  With this written explicitly, the
deduced `pi0=O(sqrt(log N/N))` and the `M` estimate are sound.  The lower
coarea proof should explicitly name the conditional second-moment
convolution above; replacing it by a change of variables that presumes
one-to-one `log R_b` would be a gap.

No claim here extends to nonregular faces, multiple receivers, or the
general strong maximal constant.

## Final read of `common_window57.md` and `geometry57.md`

Accepted, with the stated regular-window scope.  The actual full proof now
writes the needed distinction between `Pr(W<ell)` and `E g(W-ell)` in the
right order: lower concentration first makes the former polynomially small,
the target-one KKT identity identifies it with the latter, and the central
lower bound then contradicts the lower inverse-profile estimate.  From this
the displayed `pi0`, `L2` and mean-score estimates follow with the claimed
uniform powers.

The normalizer proof explicitly uses a conditional second-moment
disintegration, so it does not depend on a false global inverse for
`log R_b`.  Equations (5)--(10) have the correct signs and powers; the
`Z_+` contribution is `O(N^(3/4)L^(9/4))`, which is absorbed by the stated
`O(N^(3/4)L^3)` remainder.  The resulting strict convexity conclusion needs
only the declared connected `C^2` regular window; it correctly does not
assert existence of a stationary point or anything about `C_N`.

For the geometry note, the implicit differentiation is correct:
\(a''(0)=-K/G_{aa}\), and the fourth-order radial change is
\(K^2/(8G_{aa}G_d)\).  The local Wolfram 15.0 exact algebra was rerun and
returned `True` for this identity.  The cone-volume formula uses that `d` is
even and that the small arc has `1+d'(v)>0`; its endpoint term produces the
order-`A^4` difference, while the integral variation is only order `A^5`.
This is consistent with (11)--(12).  The saved floating diagnostic is also
internally consistent after the `Gaa` jet correction: across N=64,128,256
and both quadrature settings, `Gaa` matches `Tplus` within
`1.4e-14`, while the finite-difference Ward check is at `1.6e-9` or below.
It remains numerical corroboration only, as labeled.

## Final read of `shared57.md` — common-action hierarchy

Accepted.  For a finite action set, the loss variables lie in `[0,B]`, and
Hoeffding's lower-deviation MGF followed by a union bound gives exactly
\[
 R_J-B\sqrt{\log M/(2k)}\le \Delta_k\le R_J.
\]
The order of expectation and infimum in the argument is correct.

The three-fiber example has the displayed prefix averages.  Every pair has
a common maximizer, whereas all three have common value `8/9`; the triple
loss and the factor `2/81` in the sampled order-three certificate are
correct.  Its exact cost is `8+4 log 2`.

For the `M`-fiber construction, at prefix `j` fiber `j` alone has value
`1-delta/j`; every other fiber has value one.  A missing fiber therefore
selects a common maximising prefix for any proper selected set.  The outer
average has height `1-delta/M` only on its first cell and height one after
it, yielding `S_J=1-delta/M^2`.  The cost
\[
 M^2-\delta+(M-1)(1+\delta)\log(1+\delta)
\]
is correct: there are `M` downward cells and `M-1` upward cells, so their
linear changes total `-delta`.

Wolfram 15.0 was actually run on the staged `code/shared57.wl`.  It returned
all 34 exact checks true for `M=2,3,5,8,12` and both listed deltas, as well
as the exact three-fiber values.  This confirms the finite examples; the
arbitrary-arity assertion and Hoeffding estimate rest on the analytic
argument.  Scope remains fixed outer interval / arbitrary bounded fibers,
not a strong-maximal conclusion.
