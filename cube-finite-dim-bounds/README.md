# Fixed finite-dimensional centered cube weak-type bounds

Scope: n=3,4 first; Lebesgue background, centered axis-parallel cubes with one radius, all positive scales and general L1 inputs. Read STATE.md and rounds/005.md. Round 5 proves C3 > 1.7940007113851367; C4 lower and both upper bounds are unchanged.

Minimal independent verification of the new three-dimensional certificate (Python standard library only, no optimization and no writes):
```
python3 cube-finite-dim-bounds/audit/007c_check.py --author-result cube-finite-dim-bounds/lower/n3/round005c_certificate_result.json
python3 cube-finite-dim-bounds/src/check_round005_integration.py
```
The general L1 transfer and final geometric audit are in audit/007c_audit.md. The complete exact interval, not the shortened display decimal, is the next comparison target. MATLAB floating checks are supplied in lower/n3/round005c_matlab.m; the strict comparison uses exact Python integers/Fraction. The independent checker shares the saved input and standard integer arithmetic, but has its own root, atom-enumeration and union algorithms.

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
