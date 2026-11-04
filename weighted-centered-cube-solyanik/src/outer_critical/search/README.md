# Atomic weighted weak-type searches

This directory archives the two exploratory searches for the one-dimensional weight \(|x|^{\tau-1}\). It contains the scripts, the best phase-1 rational seed for each of 25 parameter pairs, and a compact numerical grid. It does not contain the original checkpoints or logs. Floating-point optimization is **not** a proof of optimality or a global upper bound. Exact candidate certification is handled separately.

## Dependencies and quick checks

Use Python 3 with NumPy and SciPy installed. From any working directory, with `SEARCH` set to this directory:

```sh
python3 "$SEARCH/run.py" --help
python3 "$SEARCH/phase2/run_phase2.py" --help
python3 "$SEARCH/phase2/run_phase2.py" --check-seeds
```

The seed check reads `phase1_seeds.json`, expands coalesced atoms into slots without changing the measure, and compares the Fraction-based continuous-block and all-subsets halo constructions for all 25 seeds. The scripts set BLAS/OMP thread limits to one before importing NumPy/SciPy.

## Reproduce the searches

These commands run CPU-intensive searches, with at most 15 worker processes:

```sh
python3 "$SEARCH/run.py" --workers 15 --maxiter 600 --popsize 12
python3 "$SEARCH/phase2/run_phase2.py" --workers 15 --maxiter 2000 --popsize 24
```

Each script creates its own `checkpoints/` directory beside the script and reuses existing per-task files. To start a new run, first move or remove that script's checkpoints. Phase 2 uses the 25 archived seeds in `phase1_seeds.json`; it does not require phase 1 checkpoints, and rerunning phase 1 does not automatically alter this seed file. The scripts print task failures; inspect output and checkpoints to confirm all tasks completed.

Phase 1 ran three seeds per \((\tau,N)\) for 75 tasks and performed 4,427,208 differential-evolution evaluations. Phase 2 ran 25 tasks with warm starts from a uniform chain, a two-atom configuration with small positive extra masses, and the archived phase-1 winner. Its 9,679,608 DE evaluations and 75 additional initial-seed checks bring the total to **14,106,816 DE evaluations plus 75 seed checks**. Three initial phase-2 tasks failed before optimization because phase-1 rationalization had coalesced atoms; the corrected initializer splits mass at the same position, preserving the measure, and those three tasks were retried. The other 22 checkpoints were reused.

## Scope and interpretation

The search fixes one atom at 1, constrains other positions to [−1,1], rejects positions within 10⁻⁷ of zero, and bounds positive masses by \(2^{-12}\) and \(2^5\). It enumerates continuous atom blocks to form halo intervals and uses stable power differences for the floating-point weight integral. Fraction arithmetic checks a rationalized candidate's geometry against all nonempty atom subsets; both checks use the same atomic-halo principle.

`results_summary.json` preserves the phase-2 grid, including which values are initial seeds rather than optimization improvements. In particular, the uniform 8-atom chain supplies the τ=1 value 23/16; the phase-2 maxima for τ=3/4 and 9/10 also came from the initial chain and lie below the known Melas benchmark. For τ=1/4 and 1/2 the simple symmetric two-atom configuration at ±1, with mass 4/3 each, supplies the observed best values. These are finite searches over a bounded parameter family, not evidence of global optimality. Printed numerical comparisons are not exact real-number certificates.
