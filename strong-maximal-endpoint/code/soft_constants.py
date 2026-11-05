"""Exact rational checks used to make the soft-ramp n>=1024 threshold explicit."""
from fractions import Fraction as F
import json
from pathlib import Path


def exp_bounds(x, terms=40):
    """Positive Taylor series with geometric upper bound on the remainder."""
    assert x >= 0
    power = 1
    while x > 1:
        x /= 2
        power *= 2
    term = total = F(1)
    for k in range(1, terms + 1):
        term *= x / k
        total += term
    next_term = term * x / (terms + 1)
    remainder = next_term / (1 - x / (terms + 2))
    return total**power, (total + remainder)**power


checks = {
    "log1024_gt_6_9": exp_bounds(F(69, 10))[1] < 1024,
    "log1024_lt_7": exp_bounds(F(7))[0] > 1024,
    "log6_9_lt_1_932": exp_bounds(F(483, 250))[0] > F(69, 10),
    "L_over_z2_lt_2_28": F(57, 25) * (F(69, 10) - 2 * F(483, 250)) > F(69, 10),
    "sqrt7_over32_lt_0_083": (F(83, 1000) * 32) ** 2 > 7,
    "sqrt1024_over1023_lt_1_001": F(1001, 1000) ** 2 * 1023 > 1024,
    "tail_ratio_lt_one": F(1083, 1000) * F(1001, 1000) * F(57, 25) / F(5, 2) < 1,
}
assert all(checks.values()), checks
result = {"arithmetic": "exact fractions; positive exponential Taylor remainder", "checks": checks}
output = Path(__file__).resolve().parents[1] / "rounds/r003/constant_checks.json"
output.write_text(json.dumps(result, indent=2) + "\n")
print(json.dumps(result, indent=2))
