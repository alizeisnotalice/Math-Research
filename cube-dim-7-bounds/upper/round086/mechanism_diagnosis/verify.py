#!/usr/bin/env python3
"""Bounded diagnostic arithmetic; no operator result or parameter search."""
from fractions import Fraction as F
from hashlib import sha256
from itertools import product
from math import factorial
from pathlib import Path
import json

HERE = Path(__file__).resolve().parent

def log20(y):
    z = (y-1)/(y+1)
    assert 0 < z < 1
    lo = 2*sum((z**(2*j+1)/F(2*j+1) for j in range(20)), F(0))
    return lo, lo+2*z**41/(41*(1-z*z))

alo = sum((F(-1, 2)**k/F(factorial(k)) for k in range(22)), F(0))
ahi = sum((F(-1, 2)**k/F(factorial(k)) for k in range(21)), F(0))
clo, _ = log20(1/(3-4*alo))
_, chi = log20(1/(3-4*ahi))

# Exactly eight vertices of the stated 3-node relaxation, not a geometric scan.
tau = F(1, 4)
feasible = []
for a1, a2, source in product((0, 1), repeat=3):
    if (not a1 or source <= tau) and (not a2 or source <= tau):
        feasible.append(dict(A1=a1, A2=a2, P=source, deficit=3-a1-a2-source))
minimum_deficit = min(v["deficit"] for v in feasible)
joint_energy_bound = (2-2*tau)/2
naive_sum = 2*(1-tau)

sig, chi_source, delta, gamma, c = F(9, 20), F(4, 5), F(89, 200), F(3, 2), F(9, 10)
_, nonroot_price_hi = log20(1/(1-sig))
price_budget = F(299, 500)
checks = {
    "global_critical_point_inside_4_5": F(3, 5) < alo < ahi < F(5, 8),
    "global_sqrt_coefficient_interval": F(55533956, 10**8) < clo < chi < F(55533957, 10**8),
    "mu_positive_growth_rate": log20(F(2))[1] < 1,
    "global_derivative_sign_change": log20(F(2))[0] > F(1, 2),
    "power_half_row_column_efficiency": F(1)/F(1, 2) == 2 and min(F(1), F(1)/F(1, 2)) == 1,
    "power_two_row_column_efficiency": F(1)/2 == F(1, 2) and min(F(1), F(1)/2) == F(1, 2),
    "shared_column_feasible_deficit": minimum_deficit == 1 and dict(A1=1, A2=1, P=0, deficit=1) in feasible,
    "shared_column_naive_sum_invalid": naive_sum == F(3, 2) > minimum_deficit,
    "shared_column_joint_budget_valid": joint_energy_bound == F(3, 4) <= minimum_deficit,
    "nonroot_threshold_gap": delta-chi_source*(1-sig) == F(1, 200) > 0,
    "nonroot_full_H_gap": sig-delta == F(1, 200) > 0,
    "nonroot_scalar_price_margin": nonroot_price_hi < price_budget and gamma-price_budget-c == F(1, 500),
    "sqrt_threshold_fails_same_chi_sigma": sig/2-chi_source*(1-sig) < 0,
}
assert all(checks.values()), checks
receipt = {
    "status": "PASS", "registration_sha256": sha256((HERE/"REGISTRATION.json").read_bytes()).hexdigest(),
    "scope": "13 diagnostic checks only; no new Cn theorem or exclusion of all energy constructions",
    "checks": checks,
    "exact_values": {
        "sqrt_cstar_lower": str(clo), "sqrt_cstar_upper": str(chi),
        "joint_energy_bound": str(joint_energy_bound), "naive_sum": str(naive_sum),
        "nonroot_price_upper": str(nonroot_price_hi)
    },
    "finite_graph_feasible_vertices": feasible,
    "certified_decimal_sqrt_interval": ["0.55533956", "0.55533957"],
    "three_mechanisms_only": True,
}
(HERE/"verification.json").write_text(json.dumps(receipt, indent=2)+"\n")
print(json.dumps({"status": "PASS", "diagnostic_checks": len(checks),
                  "shared_column_minimum_deficit": minimum_deficit,
                  "nonroot_threshold_and_scalar_gaps": "1/200", "no_Cn_claim": True}))
