# Positive Stieltjes pilot: substantial gain, still below1e−4 numerically

The optimal-column Stieltjes correction is analytically valid and stronger
than the normalized fixed-projection variant. The one registered numerical
pilot at11110 gives the following midpoint estimates:

| Inner level | Fixed projection columns | Optimal null-space columns |
|---:|---:|---:|
|32|8.232728542728791e−5|8.369508325092428e−5|
|64|8.232695773426486e−5|8.369391770836235e−5|

The baseline matrix midpoint value on the same outer boxes is
6.460878754203372e−5. At64 inner bins the optimal-column estimate improves
it by29.5395%, and improves the fixed-column Stieltjes estimate by1.6604%.
It is83.6939% of1e−4. The small32/64 differences are numerical convergence
diagnostics, not integration error bounds or confidence intervals.

`proof.md` derives the independently optimized feasible null-space
columns and their effective inverse coefficients. The Schur denominator
is uniformly at least11/16. The Stieltjes effective direct weight decreases
with K, so the weighted least-squares minimum K gives a valid stronger
bound. Each function and column depends only on fixed conditional
parameters; retained inverse coefficients remain1/z. The effective
matrix stays between the old harmonic and arithmetic matrices, retaining
the exact time constant and all previous full-support tail controls.

The computation reused exactly the256 outer top-two boxes. Each inner
rule used positive CDF increments on[1/4,t], and computed both Stieltjes
integrals by direct positive sums. The common low cell was represented
by its exact probability0.000195627371393305 and conditional mean
0.24002910149917345. This is an approximation; its analytic convexity
error is documented but not propagated into a strict integral enclosure.
The outer integration also remains a midpoint approximation.

The smallest sampled Schur denominator was0.6878578128594313, above the
proved11/16 guard. The largest K_opt/K_fixed was0.9999544455399922. The
sampled effective weights satisfied a_eff/a≥1.00135449 and
a_eff/m≤0.99982963. These are implementation diagnostics accompanying the
analytic inequalities, not a substitute for those proofs.

PID29404, exec session89886, exited0 after9.03seconds. The120-second
checkpoint was not reached. The final checkpoint includes source and
dependency hashes and completion status; the source supports a validated
resume from saved inner-table and orbit progress. The reproducible inner
cache is excluded from the archive whitelist. No process remains running.

No additional inner level, outer partition, full-domain grid, top-four
layer, or strict new g certificate was produced. The values above do not
upper-bound the exact formula or actual G. The inherited omitted-domain
bounds alone do not control the two numerical quadratures.

Reproduction requires the locked upper_order002 helper and result:

```sh
/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python -B cube-dim-5-bounds/upper_resolvent005/resolvent_pilot.py
```

Nothing was committed or pushed; the locked004 directory was not changed.
