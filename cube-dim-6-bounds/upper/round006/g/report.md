# A / round006: model-local g certificates, global gap retained

The finite first batch has ended. Eight continuous positive-face boxes satisfy
`G > 1/60000`, including the previous minimum region and singular/boundary
neighborhoods. This is **not** a global `g >= 1/60000` certificate and gives no
new C6 bound. The old global spectral lower bound remains the applicable
global result. The general L1/all-positive-scale operator interface remains
the primary coordinator's integration responsibility.

## Actual finite computation

The real process was Python PID **12757**, terminal execution session **23708**;
it completed normally in **59.097265916 seconds**. No MATLAB was called.
It used python-flint Arb at 384 bits, not binary64. Eight geometric points were
evaluated at epsilon = 1/200, 1/100, 1/50: 24 point/epsilon pairs, of which the
S2 point `(1,3/4,3/4,3/4,3/4,3/4)` at epsilon=1/100 is one prior certified
calibration overlap; **23 pairs are new**. Its Arb interval agrees with the
old Arb256 result. The registered limits were 8 boxes, 24 pairs and 20 minutes.

Every point uses all 21 exact fan cones and the complete 46080-element signed
permutation group, or **967680 terms** counted with exact multiplicity.
Repeated nodes and node signs are decided by integer rational arithmetic;
confluent derivatives use the existing audited formula. Identical exact
divided-difference lists, and their reflected lists, are cached algebraically.
All input hashes, output balls, node counts, and timing are in
`input_manifest.json` and `local_results.json`.

The following radii are rational and chosen by Arb comparisons in increments
of 1/100000. The box is `x1=1`, `|xi-center_i|<=r` for i=2,...,6, clipped to
the positive face. The accepted epsilon maximizes r among the three tested
values. Coordinates in this table have exact denominator 1000.

|Center numerator vector|epsilon|r|Continuous lower result|
|---|---|---|---|
|(1000,1000,893,807,467,97)|1/100|267/100000|G > 1/60000|
|(1000,990,893,807,467,97)|1/100|283/100000|G > 1/60000|
|(1000,1000,1000,750,250,0)|1/100|213/100000|G > 1/60000|
|(1000,750,750,750,750,750)|1/100|2139/100000|G > 1/60000|
|(1000,1000,1000,1000,1000,1000)|1/100|21/800|G > 1/60000|
|(1000,0,0,0,0,0)|1/50|963/20000|G > 1/60000|
|(1000,1000,0,0,0,0)|1/50|461/20000|G > 1/60000|
|(1000,900,800,600,400,200)|1/50|263/25000|G > 1/60000|

The extension uses the already audited exact inequality

    G(x) >= H_epsilon(x)
         >= exp(-r/epsilon) H_epsilon(center).

Thus every table entry is a continuous-region statement. S2 and its
intersections require no extension from their complement: H_epsilon is finite,
positive and continuous there. In particular the minimum-tie center has
`H_(1/100)` in an Arb ball around **2.17744251555e-5**, and its accepted-box
lower value is around **1.66720909882e-5**, strictly above 1/60000.

## Exact coverage gap

`check_certificate.py` reran acceptance arithmetic, rational box geometry,
the prior calibration comparison, and the new ODE algebra. It deliberately
does not claim an independent replay of the complete kernel construction.

Even after every distinct permutation of the five free face coordinates,
the sum of box volumes is exactly

    9428385958741921253 / 10^25
      = 9.428385958741921253e-7.

This is an **upper bound** for the volume covered by their union, since
overlaps may occur. Consequently the uncovered face volume is at least

    9999990571614041258078747 / 10^25.

The exact point `(1,19/20,9/10,17/20,4/5,3/4)` is outside every permuted
accepted box. The uncovered set is retained explicitly as the complement
of that finite union, not discarded by numerical sampling.

There is also a simple budget diagnosis for this specific registered batch:
with only eight boxes and its radius cap 1/10, even the most favorable
unclipped symmetric volume budget is at most
`8*120*(1/5)^5 = 192/625 < 1`. Thus this eight-box local batch could never
certify the entire face solely by this acceptance scheme. No larger uniform
grid was started. This does not exclude another adaptive or analytic method.

For cube points whose normalized direction lands in a covered box, the
already audited -5 homogeneity transfers the local certificate:
`G(rho*y)=rho^-5 G(y)>=G(y)` for 0<rho<=1, with coordinate reflection and
permutation. This does not cover directions in the stated complement.

## A valid new local estimate

`sharp_integrand_ode.md` proves by an explicit Cauchy-score square completion
that the coordinate ODE coefficient improves from 129/128 to

    4225/4224 = 1.00023674242424... .

The identity in that proof was actually checked by local SymPy, returning
exact residual zero; its constants and comparisons were evaluated in Arb384.
The tensor-product ODE certificate can therefore use
`k_epsilon=65/(8*sqrt(66)*epsilon)`. The equal-side-length improvement at a
fixed vertex lower bound is only **1.003779838658...**, so this legitimate
improvement is insufficient to resolve the coverage bottleneck.

An earlier proposed factor from applying Cauchy's derivative estimate to the
pure imaginary-axis tube was withdrawn: that bound does not control the
whole circle relative to the center density. No local box uses it. The
point computation was unchanged; only this unused coefficient was corrected
in postprocessing, with the originally executed and final source hashes both
recorded in `local_results.json`.

## Reproduction and stop status

Use the existing environment:

    cd /Users/zhengzhihao/Desktop/T/Math-Research-cube6-theory/cube-dim-6-bounds/upper/round006/g
    /Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python certify_local.py
    /Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python check_certificate.py

Necessary artifacts are this report, the finite registration, the two scripts,
the new ODE proof, input manifest, raw local results, and checked certificate.
The redundant `batch.log` terminal transcript was removed after verification.
No commits or pushes were made. The batch is stopped and awaits primary
scheduling. The next substantive advance must supply a much stronger local
estimate or a global minimization argument; more isolated point values do
not close the missing global certificate.
