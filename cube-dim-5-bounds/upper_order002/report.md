# Top-two layer improves 11110, but remains below the target

The no-tilt top-two conditional formula at11110 is rigorously between
6.08940491318362e−5 and6.92276885126884e−5. Its lower endpoint exceeds the
single-maximum layer's full upper bound4.31012868806624e−5. The new layer
therefore gives a real pointwise improvement, while its own full upper
endpoint remains below1e−4. These upper endpoints bound the corresponding
analytic lower formulas, not actual G.

All120 unordered direction pairs were grouped into11 exact invariant
classes. The full positive square for each pair was integrated: distinct
bins include both label orientations, and a same-bin weight includes both
triangles exactly once. Local SymPy verified all11 full matrix determinant
and inverse identities, as well as the diagonal probability weight.

The sole registered batch used256 geometric bins from1/4 to approximately
1.3745e6. Stable CDF and truncated moments were computed with explicit
series remainders in Arb256. The first-moment series has a same-sign
k=1,k=5 prefix; the proof and code explicitly sum that finite prefix and
apply the alternating remainder only at indices257,258,259 and beyond.
The inverse-moment multiplier has its separate conservative ratio bound.

The omitted low and high contributions are bounded above by7.226e−52
and6.763e−7. The positive full-box upper/lower bounds and exact probability
closure are detailed in `proof.md`; no numeric stable CDF or sampled
minimum is used as a certificate.

The actual run took4.48seconds, with exec session5443 terminating at
exit0. No partial checkpoint was needed. There is no running handle,
no top-three computation, no spatial grid, and no new lower-bound process.
All files remain uncommitted for parent review.

Reproduce from this independent clone:

```sh
/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python -B cube-dim-5-bounds/upper_order002/top2_certificate.py
```
