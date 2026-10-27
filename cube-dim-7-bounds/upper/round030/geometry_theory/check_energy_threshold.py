"""Exact old-energy obstruction and new union-cap checks at the frozen parameter."""
from fractions import Fraction as Q
from hashlib import sha256
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
SOURCE = HERE.parents[2] / "paper" / "constant_certificate.json"
raw = SOURCE.read_bytes()
certificate = json.loads(raw)
t = Q(certificate["t"])
A = (1 + 1 / t) ** 7
B_upper = Q(certificate["B_upper"])
B_lower = Q(certificate["B_lower"])
margin = 4 * A - B_upper
delta_old_single_box = Q(15625, 2187) / t**7
delta = Q(5**7 - 4**7, 3**7) / t**7
assert t == Q(3391749, 524288)
assert B_upper > 0 and margin > 0
assert t > 6 and B_lower > Q(54, 5)
assert A * B_upper < Q(121, 4)
assert 0 < delta < Q(1, 1000)
assert 1 + Q(3, 5) + Q(3, 5) ** 2 / 2 > Q(7, 4)
result = {
    "status": "PASS",
    "scope": "Exact B_upper < 4A and sufficient inequalities for union-cap deduction",
    "n": 7,
    "t": str(t),
    "A": str(A),
    "B_upper": str(B_upper),
    "four_A_minus_B_upper": str(margin),
    "strictly_positive_exact": margin > 0,
    "new_budget_delta": str(delta),
    "old_single_box_delta": str(delta_old_single_box),
    "new_budget_checks": {
        "t_greater_than_6": t > 6,
        "B_lower_greater_than_54_over_5": B_lower > Q(54, 5),
        "A_B_upper_less_than_121_over_4": A * B_upper < Q(121, 4),
        "delta_less_than_1_over_1000": delta < Q(1, 1000),
        "exp_series_at_3_over_5_greater_than_7_over_4":
            1 + Q(3, 5) + Q(3, 5) ** 2 / 2 > Q(7, 4),
    },
    "source_sha256": sha256(raw).hexdigest(),
    "display": {"A": float(A), "B_upper": float(B_upper), "margin": float(margin),
                "delta": float(delta)},
}
(HERE / "energy_threshold.json").write_text(
    json.dumps(result, ensure_ascii=False, indent=2) + "\n"
)
print(json.dumps({"status": result["status"], "display": result["display"]}))
