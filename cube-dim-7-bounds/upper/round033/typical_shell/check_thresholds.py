"""Exact threshold diagnostics for the all-dimensional analytic proof."""
from fractions import Fraction as Q
from math import factorial
from pathlib import Path
import json

checks = {
    "product_loss_below_two_thirds": Q(81, 127) < Q(2, 3),
    "face_count_threshold": Q(9, 128) < Q(1, 8),
    "captured_minus_corner_above_quarter": Q(5, 16) - Q(1, 100) > Q(1, 4),
    "exp_nine_halves_partial_sum_above_50": sum(Q(9, 2)**k / factorial(k) for k in range(6)) > 50,
    "exp_seven_tenths_partial_sum_above_2": sum(Q(7, 10)**k / factorial(k) for k in range(4)) > 2,
    "L_4096_x_base": Q(4096, 4) > 3 * 12,
    "n_square_loss_base": 4 * 4096**2 >= 192 * 64,
    "final_remainder_base": Q(64, 2**15) + Q(3970, 2**45) < Q(1, 256),
}
assert all(checks.values())
result = {"status": "PASS", "scope": "exact threshold base checks, no finite sampling of the all-n claim",
          "checks": checks}
Path(__file__).with_name("threshold_checks.json").write_text(json.dumps(result, indent=2) + "\n")
print(json.dumps({"status": "PASS", "checks": len(checks)}))
