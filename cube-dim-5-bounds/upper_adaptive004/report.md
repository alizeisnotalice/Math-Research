# Feasible LS correction: valid, with7.43% numerical gain in this pilot

The centered random-coefficient correction is valid and improves the
conditional inverse bound pointwise. On the one registered top-two
partition at11110, its numerical midpoint value is6.940895028998438e−5,
compared with6.460878754203372e−5 for the baseline matrix path. The gain is
7.4296%, and the resulting estimate is69.4% of1e−4. This is not a certified
interior integral or a new full-domain g bound.

`proof.md` was written before the computation. It proves VP=0 and
samplewise feasibility, the vanishing centered cross terms, the legal
dependence of all optimized coefficients on fixed conditional parameters,
the H=0 convention, the variance-upper-bound direction, and
1/m≤n_eff≤n. Retained directions always have their exact deterministic
inverse coefficient1/z. The improved matrix lies between the harmonic
and arithmetic conditional matrices, which gives a genuine pointwise
improvement while preserving the previously certified full-support tail
bound. No numerical gain is needed to justify this analytic direction.

The pilot used the same256 geometric bins as upper_order002. It evaluated
full5×5 matrices through numerical log determinants and linear solves,
without using the previous polynomial inverse evaluator. The actual
48-element stabilizer of11110 has11 orbits on the120 pairs; this exact
action was enumerated and used. The baseline midpoint value lies inside
the previous rigorous top-two interval, providing a useful consistency
check of normalization and the new matrix path.

The inherited omitted-domain upper bounds remain7.226e−52 for the low
region and6.763e−7 for the high region. These control only the omitted
domain: they do not supply an error bound for the midpoint integration
inside the partition. The source and result consistently label that
limitation. No strict new adaptive g certificate is claimed.

Actual PID28437 ran under exec session51870 and exited0 in0.94seconds.
The120-second checkpoint was not reached. The source saves the completed
orbit count, next batch, partial contributions, and source/dependency
hashes; `--resume` validates those hashes. No process remains running.
No second moment, additional partition, top-four layer, or spatial grid
was evaluated.

The analytic note `resolvent_proposal.md` records an uncomputed possible
extension: optimize over all centered one-variable corrections rather
than only the linear feature S−m. It is a proposal for a separately
registered pilot, not an additional result of this batch.

Reproduction requires the locked upper_order002 helper and result:

```sh
/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python -B cube-dim-5-bounds/upper_adaptive004/adaptive_pilot.py
```

Nothing was committed or pushed, and no locked order directory was changed.
