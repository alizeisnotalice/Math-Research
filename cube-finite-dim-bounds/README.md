# Fixed finite-dimensional centered cube weak-type bounds

Scope: n=3,4 first; Lebesgue background, centered axis-parallel cubes with one radius, all positive scales and general L1 inputs. Read STATE.md and rounds/001.md. No new bound in round 1.

From this repository root:
```
python3 cube-finite-dim-bounds/baseline/reproduce.py
python3 cube-finite-dim-bounds/baseline/calibrate.py
python3 cube-finite-dim-bounds/baseline/compare_public.py
```
These use only the Python standard library. `--save` regenerates the corresponding result summaries. `reproduce.py` checks only dimensions 3 and 4; do not use the inherited algebraic module's standalone main, which also computes dimension 2.
Run `baseline/calibrate_matlab.m` in local MATLAB for the numerical cross-check. Exact comparisons use Fraction, not MATLAB floating point.

Public lower-bound reproduction: obtain `lowerbound/check_exact2.py` and `lowerbound/c2lower_1p70272.json` at the source commit in results/public_lower_run.json, check both SHA-256 hashes, then run `python3 -I check_exact2.py c2lower_1p70272.json` in an external scratch directory. The expected exact ratio and rectangle count are recorded. No public source is silently replaced with a rounded number.

The source manifest hashes the inherited verifier and transfer proof. The inherited algebraic calculation and floor-root helper share arithmetic; they are not independent certifiers. No automated claim of a new bound is permitted without an exact comparison and the analytic transfer.
