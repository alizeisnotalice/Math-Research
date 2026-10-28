#!/usr/bin/env python3
"""Fixed s=4 eventual-theorem constants; no finite threshold or search."""
from fractions import Fraction as F
from math import factorial
from pathlib import Path
import json

HERE = Path(__file__).resolve().parent

def log12(y):
    assert y > 1
    z = (y-1)/(y+1)
    lo = 2*sum((z**(2*j+1)/F(2*j+1) for j in range(12)), F(0))
    return lo, lo+2*z**25/(25*(1-z*z))

e = F(1, 10000)
delta = F(1, 100000)
rlo, rhi = 1+delta, F(99, 25)-delta
P0 = 1-2*e
Ulo, Uhi = P0*(1-e)/(1+e), (1+e)/(P0*(1-e))
M = 2*(1+e)*Uhi
Ret = F(999, 1000)
b = 2*Ret*(1-e)*Ulo-2*e
Gamma = F(137, 100)*b/M
masslo, _ = log12(rhi/rlo)
Qlo, Qhi = log12(F(8, 5))
_, pricehi = log12(F(5, 2))
nu, eta = F(1, 4)-F(5, 4)*e, F(1, 198)
cT = F(3, 5)/2-F(7, 10)*F(2, 5)
checks = {
    "new_geometry_radius_and_ell_floor": 1 < rlo < rhi < 4 and (1+F(1, 4))/2 == F(5, 8),
    "new_Q_less_than_half": 1+F(1, 2)+F(1, 8) > F(8, 5) and Qhi < F(1, 2),
    "new_mu4_less_than_065": Qlo > F(47, 100) and 3-5*Qlo < F(13, 20),
    "new_J_derivative_constants": 2*(F(11, 5)*4)+F(3, 5)*(4+5)+1 == 24 and (24+12*4)/2 == 36 and 2*36 < 100,
    "new_J1_derivative_constants": 1+4*F(1, 2)*(F(28, 5)+F(32, 5))+4 == 29 and 29 < 40,
    "new_face_mean_rounding": F(67, 100)/(1-e) < F(671, 1000) < F(7, 10),
    "new_face_variance_rounding": 4/(3*(1-e)) < 2,
    "new_target_mass": (1-e)*masslo > F(137, 100),
    "new_weighted_cap_count_rate": nu > 0 and eta > 0 and nu*eta*F(2, 5)/(2*(1+e)) > 0,
    "new_tau_rate": cT == F(1, 50) > 0,
    "new_deficit_above_one": 0 < P0 < 1 and M > 2 and b > 0 and (1-e)*Gamma > F(13, 10) > 1,
    "new_exact_price_and_complete_slack": pricehi < F(23, 25) and F(13, 10)-F(23, 25)-F(3, 8) == F(1, 200),
    "new_full_cost_comparison_to_R080": F(3, 8) > F(28029, 100000)+F(4, 4096),
}
assert all(checks.values()), checks
values = dict(rlo=rlo, rhi=rhi, P0=P0, Ulo=Ulo, Uhi=Uhi, M=M,
              b=b, Gamma=Gamma, paid_Gamma=(1-e)*Gamma,
              mu4_upper=3-5*Qlo, target_mass_lower=(1-e)*masslo,
              log_price_upper=pricehi, nu=nu, eta=eta, tau_rate=cT)
receipt = {
    "status": "PASS", "checks": checks,
    "scope": "13 fixed-constant gates only; eventual threshold not computed; proof pending independent review",
    "exact_values": {k: str(v) for k, v in values.items()},
    "display_only": {k: float(v) for k, v in values.items()},
    "no_parameter_grid": True,
}
(HERE/"verification.json").write_text(json.dumps(receipt, indent=2)+"\n")
print(json.dumps({"status": "PASS", "fixed_gates": len(checks),
                  "paid_Gamma": float((1-e)*Gamma), "gamma": "13/10", "c": "3/8"}))
