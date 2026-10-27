"""Small exact arithmetic diagnostics, not a numerical proof or production batch."""
from fractions import Fraction as Q
from pathlib import Path
import json

epsilon0 = Q(1, 2**20)
dmax = Q(1, 2**24)
d0 = Q(1, 2**30)
variation = 1280 * epsilon0
assert variation == Q(5, 4096) and variation < Q(1, 8)
assert Q(1, 2) * Q(1, 2) * Q(1, 2) == Q(1, 8)
assert Q(1, 2) * Q(4, 3) * 2 * 4 < 8
assert Q(3, 4) * Q(1, 8) == Q(3, 32)
assert Q(1, 16) * (Q(3, 64) - Q(3, 128)) / 32 == Q(3, 65536)
assert Q(1, 65536) < Q(3, 65536)
assert 12*dmax*dmax <= dmax/Q(64)
assert 32*dmax*dmax < Q(1, 2)
delta_per_n = d0/Q(65536)
scalar_per_n = 64*d0*d0
assert delta_per_n == Q(1, 2**46)
assert scalar_per_n == Q(1, 2**54)
assert scalar_per_n < delta_per_n/Q(4)
assert delta_per_n/Q(2) == Q(1, 2**47)
result = {
    "status": "PASS",
    "scope": "exact proof-constant arithmetic only; all asymptotic/geometric statements require the analytic proof",
    "matching_log_variation": str(variation),
    "energy_deficit_per_nd": "3/65536",
    "used_deficit_per_nd": "1/65536",
    "fixed_d": str(d0),
    "operator_delta_per_n": str(delta_per_n),
    "scalar_buffer_cost_per_n": str(scalar_per_n),
    "proposed_final_gain_per_n": str(delta_per_n/2),
}
Path(__file__).with_name("constant_checks.json").write_text(json.dumps(result, indent=2) + "\n")
print(json.dumps({"status": "PASS", "scope": result["scope"]}))
