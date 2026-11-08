# Final audit of `tail71.md` and `cell_certificate71.md`

## Verdict

The two notes are mathematically coherent in their stated scopes. The
fixed-`N` tail, the noncommuting `D_N` limits, the compact-column
`O(N^{-1})` refinement, and the certified-cell construction are valid in
their displayed parameter ranges. None supplies a
dimension-uniform weak `(1,1)` estimate.

## `tail71.md`

### Fixed-dimensional tail: PASS, with scope made explicit

For nonnegative `G` of finite physical mass supported in the stated box, the
tail proof has the right direction. At a source--receiver separation with
`T=-\sum x_i`, the full counted-label sum gives

\[
 H_{N,c}(T+w)\sim {2^N\over\Gamma(N)}T^{N-1}.
\]

On the outer cone all increment labels are allowed, so their permutation of
the source labels yields exactly the counted mass `m`, not an extra `2^N`.
The counted receiver-label shell has the compensating `2^N`. Thus the tail
constant `m` is correct under the counted convention; under Haar
normalization both factors disappear.

The `q\pm O(q^{-1})` cutoff is also justified. A relative field error
`1+O(q^{-1})` is compensated by changing `T` by `C/q`, since

\[
 -\partial_T\log a(T)=1+(N-1)/T\ge1.
\]

Choose `C` larger than the uniform coefficient in the field error. Then
the lower outer cumulative cone with `T\le q-C/q` is inside the strict
superlevel set, while that set is contained in the reachable upper
cumulative cone with `T\le q+C/q`. The inclusions have the correct
direction.

The cone-volume comparison has the right orientation. For the reachable
upper cone `x_i\le R`, its polynomial is `(T+NR)^{N-1}`; for the outer
lower cone `x_i\le-R-c`, it is `(T-N(R+c))^{N-1}`. Replacing the first
`NR` by `N(R+c)` is a harmless fixed-`N` `1+O(T^{-1})` enlargement. In the
display after (3), the fixed cone shifts are absorbed into the already
written `1+O_{N,c,R}(q^{-1})`; this is sufficient for the limiting
conclusion.

The first section needs the explicit hypothesis `G_\sigma\ge0` (as is
natural for the local-cost use). Mere finiteness of the signed quantity `m`
would not give either comparison in (2).

### Critical column: PASS

For fixed `N,c`, eventual polynomiality gives `D_N(s)\to1`. In the
critical sequence, `\partial_s\log H=1+o(1)` yields the fixed-`r` ratio
`H(s-r)/H(s)\to e^{-r}`; monotonicity gives the integrable domination
`e^{-r}` after the outside factor in `D_N`. Hence the limit is `1/2`, and
the two limiting procedures genuinely differ.

The compact refinement is also sound. From `B=\gamma+O(N^{-1})` and
`\gamma'=O(N^{-1})`, without differentiating the remainder, one gets

\[
 \int_0^tB(s-v)dv=B(s)t+O((t+t^2)/N),\qquad0\le t\le\sqrt N.
\]

Although this exponent error is only `O(1)` at the endpoint, it is bounded
there, and integration against `e^{-t}` makes its contribution `O(N^{-1})`;
the `t>\sqrt N` tail is exponentially small. The expansion at `c=\mu`
then has the displayed sign and coefficient. Its stated conditional,
macroscopic scope must remain in force.

### Posterior truncation: PASS

The log-concavity ratio bound and the Gaussian-in-index tail are correct.
Keeping all label configurations with counts `|K-m|\le R_{\rm cut}` and
renormalizing changes any `[0,1]` test by at most the omitted posterior
mass, so (15) is valid (and applies to nonsymmetric tests as well, provided
the conditional distribution within each retained count is preserved).

The current internal references are correct: (14) follows by iterating
(13), and (15) uses (14). The binomial factor represents all
`\binom Nk` label configurations of each count; the conclusion for arbitrary
`[0,1]` tests is valid provided the conditional label distribution within
each retained count is preserved. The lemma implicitly needs `s>0`, as is
automatic for the conditional-sum setting.

## `cell_certificate71.md`: PASS with its stated label restriction

The Jacobian identity is correct: `dz=ds\,dH/\sqrt N` turns the conditional
factor `e^{-s}/(\sqrt N H)` into `e^{-S}/H` in Cartesian increments. The
decreasing base density is therefore minimized at the lower receiver corner
`v=j\Delta`, producing the exact Irwin--Hall coefficient `q_M`.

The two certificate restrictions are necessary and sufficient for the
displayed retained lower contribution:

\[
 k_i\ge j_i+1,\qquad
 c=d\Delta\implies J_i=1\text{ uniformly only if }k_i-j_i\ge d+1.
\]

Thus `2^{\#\{i:k_i-j_i\ge d+1\}}q_M` is sound. A floating coefficient
using `k_i-j_i\ge d`, or allowing zero coordinate gaps, cannot certify a
cellwise lower bound.

One important scope condition is already present and should not be weakened:
the label multiplier is pointwise valid for label-invariant source cells (or
for a common lower height over every source label). It is not a lower bound
for a general source merely from its label-average. This is harmless for a
weak-type *lower-bound construction*, since one may choose label-invariant
test inputs.

Finally, the interval directions and strict margin in (10) are correct:
use lower matrix/receiver-mass enclosures, upper source-mass enclosures, and
a positive margin above `\lambda`. They turn selected disjoint half-open
cells into a genuine strict-superlevel witness. No current r070 floating
MILP output meets those certification requirements.

