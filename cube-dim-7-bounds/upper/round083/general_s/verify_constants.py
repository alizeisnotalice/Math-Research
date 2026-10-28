#!/usr/bin/env python3
"""Independent fixed s=9/2 constants; no finite threshold or producer replay."""
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import argparse
import json

HERE = Path(__file__).resolve().parent
parser = argparse.ArgumentParser()
parser.add_argument("--constants-dir", type=Path, default=HERE.parent/"constants")
args = parser.parse_args()
reg = json.loads((HERE/"REGISTRATION.json").read_text())
for name, digest in reg["producer_locks"].items():
    assert sha256((args.constants_dir/name).read_bytes()).hexdigest() == digest
producer = json.loads((args.constants_dir/"result.json").read_text())
P = {k: F(v) for k, v in reg["tuple"].items()}

def log12(y):
    assert y > 1
    z = (y-1)/(y+1)
    lo = 2*sum((z**(2*j+1)/F(2*j+1) for j in range(12)), F(0))
    return lo, lo+2*z**25/(25*(1-z*z))

s, e, qmax = P["s"], P["e"], P["qmax"]
rlo, rhi = 1+P["delta"], s*qmax-P["delta"]
Qlo, Qhi = log12(F(18, 11))
muhi = F(7, 2)-F(11, 2)*Qlo
masslo, _ = log12(rhi/rlo)
_, pricehi = log12(1/(1-P["sigma"]))
old81price, _ = log12(F(5, 2))
P0 = 1-2*e
Ulo, Uhi = P0*(1-e)/(1+e), (1+e)/(P0*(1-e))
M = 2*(1+e)*Uhi
b = 2*P["Ret"]*(1-e)*Ulo-2*e
Gamma = P["mL"]*b/M
nu = P["qmin"]-(1+P["qmin"])*e
eta = (1-qmax)/(2*qmax)
tau_rate = P["sigma"]/2-P["chi"]*(1-P["sigma"])
A_prime = 2*F(25, 11)*F(9, 2)+F(7, 11)*(F(9, 2)+F(11, 2))+1
J_prime = (A_prime+F(27, 2)*F(9, 2))/2
J1_prime = F(9, 8)+F(9, 4)*(F(72, 11)+F(81, 11))+F(81, 16)
target_bad_constant = F(6561, 256)*50000**2
checks = {
    "producer_identity_parsed_not_replayed": producer["status"] == "PASS" and len(producer["gates"]) == 18 and all(producer["gates"].values()),
    "new_geometry_domain": s == F(9, 2) and P["qmin"] == 1/s and 0 < P["qmin"] < qmax < 1 and 1 < rlo < rhi < s,
    "new_relative_gap_and_Q_bound": (1+1/s)/2 == F(11, 18) and F(79, 48) > F(18, 11) and Qhi < F(1, 2),
    "new_J_derivative_arithmetic": A_prime == F(306, 11) and J_prime == F(3897, 88) and J_prime < 45 and 2*45 == 90,
    "new_J1_derivative_arithmetic": J1_prime == F(6597, 176) and J1_prime < 40 and s/4 == F(9, 8),
    "new_J2_and_baseline_constants": s/F(12) == F(3, 8) and F(27, 4) < 7 and 7**2 == 49 and 4*7**2 == 196 and 2*7 == 14,
    "new_mu_bound": 0 < muhi < F(7914, 10000) < P["target"],
    "new_target_probability_constant": target_bad_constant == 64072265625 and F(9, 8)*s*100000 == 506250,
    "new_transport_gap": P["transported"]-P["target"] == F(1, 100000) and 2*40 == 80 and 2*90 == 180,
    "new_inverse_face_mean_and_gap": P["transported"]/(1-e) < P["face"] and P["chi"]-P["face"] == F(3, 100000) > 0 and Qhi < P["pin"],
    "new_true_variance_bound": F(27, 16)/(1-e) < 2,
    "new_weighted_cap_count": nu == F(22221, 100000) > 0 and eta == F(1, 19998) > 0,
    "new_tau_rate": tau_rate == F(379, 2000000) > 0,
    "new_target_mass": (1-e)*masslo > P["mL"],
    "new_deficit_above_one": 0 < P0 < 1 and M > 2 and b > 0 and (1-e)*Gamma > P["gamma"] > 1,
    "new_price_and_complete_slack": pricehi < P["price"] and P["gamma"]-P["price"]-P["c"] == P["slack"] == F(1, 1000),
    "old81_exact_positive_price": old81price > F(229, 250),
    "full_cost_comparison_old81": P["c"] > F(13, 10)-F(229, 250)+F(5, 64),
    "full_cost_comparison_old80": P["c"] > F(39029, 100000)-F(11, 100)+F(4, 64),
}
assert all(checks.values()), checks
values = dict(rlo=rlo, rhi=rhi, Q_upper=Qhi, mu_upper=muhi,
              J_prime_bound=J_prime, J1_prime_bound=J1_prime,
              target_bad_constant=target_bad_constant,
              target_mass_lower=(1-e)*masslo, P0=P0, M=M,
              Gamma=Gamma, paid_Gamma=(1-e)*Gamma,
              face_mean_after_inverse=P["transported"]/(1-e),
              face_gap_rate=P["chi"]-P["face"],
              weighted_cap_nu=nu, weighted_cap_eta=eta,
              tau_rate=tau_rate, price_upper=pricehi,
              old81_price_lower=old81price,
              old81_comparison_margin=P["c"]-(F(13, 10)-F(229, 250)+F(5, 64)),
              old80_comparison_margin=P["c"]-(F(39029, 100000)-F(11, 100)+F(4, 64)))
receipt = {
    "status": "PASS", "registration_sha256": sha256((HERE/"REGISTRATION.json").read_bytes()).hexdigest(),
    "recipe": reg["independent_recipe"], "checks": checks,
    "exact_values": {k: str(v) for k, v in values.items()},
    "display_only": {k: float(v) for k, v in values.items()},
    "scope": "19 new fixed-constant gates; no threshold or operator acceptance; whole-chain proof separate",
}
(HERE/"verification.json").write_text(json.dumps(receipt, indent=2)+"\n")
print(json.dumps({"status": "PASS", "new_fixed_gates": len(checks),
                  "paid_Gamma": float((1-e)*Gamma), "gamma": str(P["gamma"]), "c": str(P["c"])}))
