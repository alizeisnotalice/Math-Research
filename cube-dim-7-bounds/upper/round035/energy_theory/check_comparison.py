"""Exact comparison of two frozen budget coefficients; no search or optimization."""
from fractions import Fraction as Q
from pathlib import Path
import json

q = Q(74275164113, 10**12)
threshold = (65536 * q)**2
below = threshold.numerator // threshold.denominator
above = below + 1
assert Q(below) < threshold < Q(above)
assert Q(below, 65536**2) < q*q < Q(above, 65536**2)
result = {
    "status": "PASS",
    "scope": "exact comparison only; does not supply the R035 asymptotic theorem threshold",
    "old_q": str(q),
    "log_log_n_strict_comparison_threshold": str(threshold),
    "integer_below_threshold": below,
    "integer_above_threshold": above,
    "both_theorems_must_also_be_applicable": True,
}
Path(__file__).with_name("comparison_check.json").write_text(json.dumps(result, indent=2) + "\n")
print(json.dumps({"status": "PASS", "integer_bracket": [below, above]}))
