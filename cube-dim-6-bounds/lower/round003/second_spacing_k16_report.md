# Round 003 local-line winner: primary exact certificate

The pre-registered line `a_2 -> (1+k/10000)a_2`, `k=2,...,16`, completed
with no background process.  It recomputed the finite witness union at all
15 points; all had 127 raw and 96 noncontained boxes.  The largest floating
diagnostic was `k=16`, approximately `2.0656325551638566`.  The diagnostic
only selected the candidate.

The single permitted promotion was then run exactly at `k=16`.  In
`Q(u,v)`, `u^6=2,v^6=3`, it regenerated the boxes, checked each coordinate
endpoint against both block extremes at the common radius, removed
contained boxes, and recursively computed the full six-dimensional union.
The resulting exact interval has rational enclosure

```
206563255516385/10^14 < R_(k=16) < 206563255516386/10^14.
```

The calculation also made the exact comparison that its lower endpoint
strictly exceeds the full frozen R6-001 upper endpoint.  The primary
certificate is [pending an independent non-importing audit]; it imported
the round-001 number-field layer only for arithmetic, so it alone is not a
new certified bound.  See `certify_second_spacing_k16.py` and
`second_spacing_k16_certificate_results.json` for the complete interval,
all 96 witness count vectors, and the 36 exact coefficients.

The finite-L1 transfer proof in `second_spacing_transfer.md` applies after
replacing `10001/10000` by `10016/10000`: it uses compact inner witness
boxes, small equal-mass smooth bumps, source collar two lattice cells, and
the limits in block size, smoothing width, and strict threshold.  This
keeps the fixed Lebesgue/common-radius operator unchanged.
