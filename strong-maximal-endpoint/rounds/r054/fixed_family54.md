# r054 fixed-family certified comparison

## Inputs and scope

This archive fixes the r052 geometry `kappa=-0.08`, the 65 mirror-paired
receivers, the cap `2`, target `1`, and the finite Pareto-cone witness.  The
immutable interval core has SHA-256
`17ef3f7fbd83211f7a411a6700ba751a94a1d5a85bb1c3b7dff9333d968bc3de`.
The D1 record and cone record capture their core, r052 input, r030 seed, and
final re-read hashes; equality of startup and final hashes is required.
The historical absolute `source_snapshot` pathname inside the D1 replay JSON
is not a portability requirement: this archive's `cert_d1_core54.py` is
byte-for-byte equal to that source, as recorded in `manifest54.json`.

No statement here is about the full strong-maximal halo, arbitrary receiver
subsets, arbitrary geometries, or the constant `C_N`.

## D1 lower bound

The frozen core encloses the finite Gamma--Bernoulli mixture and a concave
chord lower minorant of

\[
q(a)=a-\Phi_{64,[0,2]}^*(a).
\]

The window contribution is evaluated from exact CDF increments and the
outside-window loss uses `q(a)>=-a` and a union bound.  The replay record
therefore certifies

\[
D_1\ge d_-=
0.0551647563805416712476194837891016143171282088410593889229509027188992503482270807998388018.
\]

The displayed narrow interval called `minorant_window_functional_interval`
in the JSON is an interval for the **chosen chord minorant functional**.  It
is not an upper bound for the actual D1 functional.

## Finite cone volume

For the 65 fixed anchors, order the distinct first coordinates as
`x_1<...<x_m` and put
`b_i=min{y_j:x_j<=x_i}`.  The finite Pareto-cone union has the exact
decomposition

\[
V=\sum_{i=1}^{m-1}(e^{-x_i}-e^{-x_{i+1}})e^{-b_i}
 +e^{-x_m}e^{-b_m}.
\]

All exponents are exact finite decimals.  The frozen core's outward `exp`
enclosures give

\[
 V\le 17.3162397952482031986842475186818696906396407065561880718529347521923696079925953106117130.
\]

## Weak-dual implication

Let `H` satisfy `0<=H<=2` and every one of the fixed 65 normalized receiver
averages be at least one.  The capped Fenchel inequality says

\[
\Phi_{64}(H)\ge aH-\Phi_{64,[0,2]}^*(a).
\]

For this fixed mixture, the 65 average constraints imply the linear term
integrated against the mixture is at least its target-one value.  Integrating
the inequality consequently gives

\[
J\ge e^{-\tau}d_-,\qquad \tau=65.68292727325189.
\]

Hence this one fixed witness has score at most

\[
\tau+\log V-\log d_-
\le 71.4320030519008984021498728220331218985271992312872065540724419851339496073718634077432033
<\frac{143}{2}=71.5.
\]

This only applies under the stated boundedness and all-65-average hypotheses.
It does **not** prove that a nonmonotone `H` generates these cones, control
extra halo points, or extend to any receiver subfamily.  Any comparison with
an independently certified scalar lower bound must be made by that separate
certificate; this archive supplies only the fixed-family upper side.

## Interpolated old-anchor Pareto chain

The exact 65 anchors, sorted by their first coordinate, were checked to have
strictly increasing first coordinate and strictly decreasing second
coordinate.  Join consecutive old anchors by straight segments.  With
`dx=x_1-x_0>0`, `dy=y_1-y_0<0`, and `z=dx+dy`, direct parameterization gives

\[
 \int_0^1 dx\,e^{-(x_0+t\,dx)-(y_0+t\,dy)}\,dt
 =dx\,e^{-x_0-y_0}F(z),\qquad
 F(z)=\frac{1-e^{-z}}z,\quad F(0)=1.
\]

Adding the last tail `exp(-x_last-y_last)` gives the continuous-cone volume

\[
 V_\Gamma\le
21.9276177716117822387060492721956642698526306379295403308967416744335070327342282793885916.
\]

The interval code handles either sign of `z` by a general directed interval
division; in this chain the exact range was
`[-0.0864020752727590,0.0864020752727590]`.  Reusing the same old-65 D1
lower bound gives

\[
 \tau+\log V_\Gamma-\log d_-
 \le71.6681052020114122122248411656321424198417538368233986472272008940177404560085488911287349
 <71.79.
\]

This is an exclusion statement only for a proposed refined receiver family
whose witness is this particular interpolated chain **and** whose input
still satisfies all old 65 target-one constraints.  New receiver parameters
along the chain may be arbitrary, but cannot replace the old constraints in
this implication.  It neither constructs a halo for a nonmonotone input nor
controls any additional halo points.

## Reproduction

Run the portable launcher from a clone whose root contains `rounds/r052` and
`rounds/r030`:

```sh
python3 certify_d1_54.py --project-root /path/to/strong-maximal-endpoint \
  --m 512 --workers 8 --out d1_m512_replay.json
python3 cone_score54.py --project-root /path/to/strong-maximal-endpoint \
  --d1-json d1_m512_replay.json --out cone_score54.json
```

The `--m 16` launcher run is a real path/I/O/hash smoke test but deliberately
does not pass the positivity threshold; it should exit nonzero after writing
its diagnostic JSON.
