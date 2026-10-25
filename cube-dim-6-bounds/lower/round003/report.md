# Round 003 lower: rational spacing perturbation screen

The frozen comparator was the complete exact upper endpoint of `R6-001`
from the main project's current ledger, read with SHA-256
`ab862aa4711d74665d583414c28c4e304bb42484569556758756eb3c77872088`.
It was not replaced by the shorter decimal display.

All 16 pre-registered rational perturbations completed in 0.598 seconds.
Each recomputed the nonempty and noncontained common-radius witness boxes;
each had 127 and 96 respectively.  The best discovery value was candidate
4, with spacing increment
`(-1,+1,+1,+1,+1,+1)/500`, and its full-union floating diagnostic was
`2.0629276496319466`.  It is lower than the old certified display upper
by `0.0022082931337434886`, far beyond the `10^-7` promotion margin.

Consequently none was sent to exact number-field or Fraction union
certification.  This rules out only these sixteen listed diagnostics under
the fixed `ell in {1,2,3,4}^6` witness screen.  It is not an upper bound
on any perturbation family and does not refute satellite, larger-block, or
other six-dimensional constructions.

The actual command was

```
/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python -u \
  cube-dim-6-bounds/lower/round003/screen_spacing_perturbations.py \
  --ledger /Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-6/cube-dim-6-bounds/results/bounds_ledger.json
```

It exited zero and created `screen_results.json`.  No persistent process
or MATLAB session was started.
