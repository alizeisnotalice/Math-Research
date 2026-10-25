# Top-three pilot is not presently competitive with128e−6

The rigorous top-three conditional interface is valid. At11110, the one
registered numerical pilot estimates its formula at7.816013872981878e−5,
about61.1% of the target128e−6. Eight consecutive block means range from
7.80906e−5 to7.82777e−5. These are numerical diagnostics only: there is no
strict integration error or tail certificate, and the block extrema are
not bounds or confidence intervals.

The exact geometry has560 unordered triples and28 coefficient classes.
All retained Gram matrices have rank3; their determinants are48,64,112.
The determinant/inverse formulas use nonnegative Gram minors and remain
well-defined when augmented minors vanish. `proof.md` gives the event
partition, all repeated-bin weights, the13-variable truncation, and the
precise560×6 normalization for the ordered sampling representation.

The sole pilot used16384 scrambled Sobol points and a fixed4096-interval
CDF/moment table. Stable table values came from the exact-series helper
in `../upper_order002/top2_certificate.py`; its hash is recorded in the
result. Interpolation and integration were numerical, with no table
extrapolation. Both the code and dependency hashes are recorded.

Actual PID27535 ran under exec session14365 and exited0 in1.60seconds.
The120-second checkpoint was not reached. The code saves table progress
and the accumulated geometry sum, with source-hash validation for
`--resume`; the final checkpoint records completion. Checkpoint arrays
are local reproducible cache, outside the archive whitelist.

No further box partition, strict integral certification, top-two
refinement, or full-domain grid was started. This is not a proved top-three
obstruction, and it does not upper-bound actual G. The low numerical value
supports moving the next finite test toward a different loss mechanism.

`adaptive_ls_proposal.md` derives a concrete alternative: random
least-squares coefficients that remain feasible for every realization.
Independent centered corrections make the optimization diagonal and
improve the harmonic inverse coefficient, using only existing truncated
first/inverse moments and the bound Var(S)≤m(t−m). This proposal has not
yet been numerically evaluated.

Reproduce, keeping the previously locked upper_order002 helper available:

```sh
/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python -B cube-dim-5-bounds/upper_order003/top3_pilot.py
```

Nothing was committed or pushed. There are no running handles.
