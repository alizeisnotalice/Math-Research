# lower004: two single-spacing endpoint screens

## Scope and batch

This independent diagnostic starts from the current R6-003 breakpoint
construction and its complete exact upper endpoint, read from the current
ledger at runtime.  It holds the new second spacing fixed and varies only
`a1`, then only one representative tail spacing `a3`; the remaining tail
coordinates remain fixed.  It does not overlap the separate period-two mass
search.

For each axis the program formed every pre-branch endpoint function from
count vectors with nonempty other five coordinate slabs, merged positive
pair-equality values numerically for selection only, and evaluated at most
eight nearest events below plus eight above the base.  It then re-enumerated
all `4^6` count vectors and computed the full union at each selected value.
The exact event descriptors (rational slope, radius sign, mass and count) are
in the generated JSON, allowing a later algebraic checker to reconstruct a
candidate without trusting its displayed decimal.

The first screen implementation inadvertently treated a symbolic base event
as slightly below the base because two mpmath expression trees rounded
differently.  It was rejected as a repeated baseline, a symmetric `10^-50`
relative exclusion was added, and the finite batch was rerun.  The reported
result is the corrected run only.

## Actual corrected run

Executed with Python 3.15 and mpmath 1.3.0:

```
/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python \
  cube-dim-6-bounds/audit/lower004/endpoint_event_screen.py
```

It completed successfully in about 22 seconds.  The `a1` event generator had
210 distinct endpoint forms and 1,331 positive merged events before choosing
16; the tail generator had 165 forms and 883 events before choosing 16.

No one of the 32 numerical event values exceeded the complete R6-003 upper
endpoint.  The best screened values were:

| varied coordinate | best event ratio | deficit from current complete upper |
| --- | ---: | ---: |
| `a1` | `2.066553479024977523...` | `-3.703469859...e-4` |
| `a3` | `2.066913981003423294...` | `-9.845007457...e-6` |

The raw/maximal counts varied (`122/96`, `127/95`, `127/96`, `132/100` on
the `a1` batch and `126/95`, `127/96` on the tail batch).  This is recorded
as diagnostic evidence of event changes only.  Equal counts were never
used to infer the same union decomposition, and unequal counts are not by
themselves a proof of a new mechanism.

## Status

No new lower bound, exact candidate, or optimality claim results from this
round.  The registered 16-event limit per axis has been reached; this axis
screen is stopped.  The code uses the earlier Fubini/one-spacing piecewise
affine interface only to define a finite endpoint search, not to infer an
unverified interval formula at a selected event.
