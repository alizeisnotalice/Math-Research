# r061 finite-multiplier dual certificate

`certificate61.py` is a replayable directed-Decimal verifier for the fixed
33 anchors

\[
(0,0),\qquad (v+v^2/320,-v+v^2/320),\qquad
(-v+v^2/320,v+v^2/320),\quad 1\le v\le16.
\]

It uses the archived exact decimal `c=1.2360667450799148`, cap `2`, and the
rounded nonnegative finite multipliers recorded in `certificate61.py` from
the selected r060 gamma=`1/5` seed.  Its target is deliberately weaker than
the floating r060 diagnostic: it proves `e^tau P >= 0.1000` for the finite
target-one constraint system, if the emitted JSON reports `target_met: true`.

For finite multipliers \(\lambda_j=e^{-\tau}w_j\), write
\(S=\sum_j\lambda_jR_j\) and \(q(a)=a-\Phi^*_{[0,2]}(a)\).  Fenchel weak
duality and normalization of every kernel give

\[
 e^\tau P\ge e^\tau\int q(S)\,d\mu.
\]

On a rectangle in the square `[-15.2,160]^2`, `S` is coordinatewise
increasing.  A certified concave chord of `q` over the resulting interval is
therefore a lower bound.  Its affine integral is computed from exact
Gamma(32)--Bernoulli CDF increments, not numerical quadrature.  Outside the
square the proof uses `q(S)>=-S` and the union bound on the two coordinate
tails.  The emitted `tail_loss_upper` is exactly that loss.

The cone-volume field is a separate finite-anchor calculation.  Neither it
nor this finite dual certificate establishes a general endpoint-constant
claim.  It does, however, immediately apply to the stated *fixed continuous
parabola family*: every input satisfying its target-one constraints for all
`v in [-16,16]` satisfies the 33 constraints used here.  Hence its source
cost is at least the emitted finite-system lower bound (in particular at
least `0.1000 exp(-tau)`).

For that same fixed curve, its union of upper orthants has the exact value

\[
 V=e^{-8/5}+32\sum_{k\ge0}\frac{(-8/5)^k}{k!(2k+1)}.
\]

Indeed, parameterizing the decreasing boundary by `v` gives
`x=v+v^2/320`, `y=-v+v^2/320`, and therefore
\(V=e^{-8/5}+\int_{-16}^{16}e^{-v^2/160}\,dv\).  An alternating-series
upper bound through the positive degree-12 term, together with
\(e^{-8/5}<1/\sum_{k=0}^{16}(8/5)^k/k!\), proves the convenient rational
bound \(V<21\).  Consequently the certified `0.1000` lower bound yields
the concrete fixed-family score bound

\[
 \log(V/P)<\tau+\log 210<71.79.
\]

The final strict inequality can be checked entirely with rationals by the
positive degree-30 Taylor lower sum for
\(e^{71.79-\tau}\), which exceeds `210`.  This consequence remains only
about the specified continuous parabola and capped source problem; it is not
a claim for arbitrary receivers, an all-family `C`, or a strong-maximal
endpoint constant.

## Why the stored chords are lower bounds

The capped conjugate gives
`q(a)=min_{0<=h<=2}{Phi_64(h)+(1-h)a}`, so `q` is concave,
`q(a)=a` on `[0,1]`, `q=1` on `[1,64]`, and it is nonincreasing above64.
For a rational height `1<h<2`, put `a*=Phi_64'(h)`. Then
`q(a*)=Phi_64(h)+(1-h)a*`. The code chooses `x<=a*` as the directed
lower derivative endpoint and
`y=Phi_lower(h)+(1-h) a*_upper`. Since `x>=64` and `1-h<0`,
`y<=q(a*)<=q(x)`. The cap node is placed above the true derivative and
uses the exact affine cap branch, with a lower `Phi_64(2)` endpoint.
Thus every stored node is below the true `q`; chords between such nodes
are lower bounds by concavity. For a cell with certified `L<=S<=U`, the
chord between the returned lower values at `L,U` is again below `q(S)`.
Its affine integral uses directed arithmetic with sign-aware interval
multiplication, not a Riemann point sample.

Each shifted kernel is coordinatewise nondecreasing as a likelihood
relative to `exp(-t1-t2)dt`. The endpoint polynomial sums therefore bound
`S` throughout a rectangle. Its affine integral's `S` term is exactly a
weighted sum of products of one-dimensional Gamma--Bernoulli CDF
increments. All input endpoints and weights are exact finite decimals;
`exp` and `ln` are enclosed by correctly rounded Decimal evaluations and
one outward adjacent representable endpoint. No optimality claim about
the chosen multiplier is required.

Every kernel has support above its own anchor, hence there is no omitted
lower-coordinate probability below `LO=-15.2`. The upper-coordinate union
bound controls the remaining tail; `|q(S)|<=S` ensures all integrals exist.
The use of reflection symmetry in the rectangular sum is exact because
both the grid and the weighted anchor list are invariant under swapping
the two coordinates.

Run from the project root:

```sh
python3 code/certificate61.py --m 256 --out /tmp/certificate61.json
```

The completed run gives a directed lower bound greater than `0.1128=141/1250`. Combining it with the rational Taylor bounds in geometry61 gives the fixed full parabola-cone score below `71.666`. This is a one-way consequence for inputs satisfying every continuous-curve constraint, not a certificate that an input meets them. The displayed upper interval endpoints enclose only the selected chord-minorant integral, not the unknown primal cost. Any nonnegative finite rational multipliers are valid; proximity to the r060 optimizer is not a hypothesis.
