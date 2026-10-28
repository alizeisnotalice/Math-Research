#!/usr/bin/env python3
"""Independent R089 compressed threshold audit; no producer guard execution."""
from fractions import Fraction as F
from hashlib import sha256
from math import factorial
from pathlib import Path
import argparse
import json

HERE = Path(__file__).resolve().parent
parser = argparse.ArgumentParser()
parser.add_argument("--source-dir", type=Path, required=True)
parser.add_argument("--generic-proof", type=Path, required=True)
parser.add_argument("--main-acceptance", type=Path, required=True)
parser.add_argument("--full-review", type=Path, required=True)
args = parser.parse_args()
reg = json.loads((HERE/"REGISTRATION.json").read_text())
for fn, digest in reg["producer_locks"].items():
    assert sha256((args.source_dir/fn).read_bytes()).hexdigest() == digest
assert sha256(args.generic_proof.read_bytes()).hexdigest() == reg["generic_R088_proof_sha256"]
dependency = json.loads((HERE/"DEPENDENCY_ACCEPTANCE.json").read_text())
assert sha256(args.main_acceptance.read_bytes()).hexdigest() == dependency["main_acceptance_sha256"]
assert sha256(args.full_review.read_bytes()).hexdigest() == dependency["independent_full_review_sha256"]
passed = json.loads((args.source_dir/"candidate2_result.json").read_text())
failed = json.loads((args.source_dir/"candidate_result.json").read_text())

def logs8(y):
    assert y > 1
    z = (y-1)/(y+1)
    lo = 2*sum((z**(2*j+1)/F(2*j+1) for j in range(8)), F(0))
    return lo, lo+2*z**17/(17*(1-z*z))

def T8(y):
    assert y > 0
    return sum((y**j/F(factorial(j)) for j in range(9)), F(0))

# Exact parameters independently transcribed from the locked fixed registration.
L = F(2**27)
s, sig, chi, dt = F(5), F(5013, 10000), F(201, 200), F(12532, 25000)
qmax, trim, e = F(49999, 50000), F(1, 10000000), F(1, 100000)
Ret, mass, gamma, face = F(999, 1000), F(4007, 2500), F(3193, 2000), F(247, 250)
price, c, allowance, pin = F(3479, 5000), F(9, 10), F(7, 10000), F(3, 5)
x, h, N = (1-sig)*L, 1/L**2, L**12/F(factorial(12))
rlo, rhi, qmin = 1+trim, s*qmax-trim, 1/s
rho = N/(s*x)
Qlo, Qhi = logs8(F(5, 3))
mulo, muhi = s-1-(s+1)*Qhi, s-1-(s+1)*Qlo
dlo, dhi = (chi-muhi)/4, (chi-mulo)/4
masslo, _ = logs8(rhi/rlo)
_, pricehi = logs8(1/(1-sig))
P0 = 1-2*e
Ulo, Uhi = P0*(1-e)/(1+e), (1+e)/(P0*(1-e))
M = 2*(1+e)*Uhi
b = 2*Ret*(1-e)*Ulo-2*e
Gamma = mass*b/M
mean = 81*s*s*x*x/N+18*s*x/N
variance = 324*s*s*x*x/N
nu, eta = 1/s-(1+1/s)*e, (1-qmax)/(2*qmax)
K = (nu-1/N)*eta*x/(1+e)
face_gap = (chi-face)*x-pin
face_fail = 10*x/face_gap**2
oldx = (1-sig)*F(2**26)
old_face_fail = 10*oldx/((chi-face)*oldx-pin)**2
ct = dt-chi*(1-sig)
cap = M/T8(K)
tau = rhi*L**2/T8(ct*L)
scalar = 1/N+1/T8(dt*L)+1/T8(sig*L)+L/T8((sig-dt)*L)
target_cap_fail = F(19, 20)**511
checks = {
    "source_identity_and_historical_status": passed["status"] == "PASS" and len(passed["gates"]) == 56 and all(passed["gates"].values()) and failed["status"] == "FAIL" and failed["failed"] == ["face_retention"],
    "fixed_geometry_domain": s == 5 and qmin == F(1, 5) and qmin < qmax < 1 and 0 < sig < 1 and 0 < e < F(1, 12) and N >= 512 and 1 < rlo < rhi < 5,
    "complete_active_matching": 2*h < trim and rlo-2*h > 1 and 1-h > 0 and (s-1)*qmin-h > 0 and h < F(1, 1000),
    "finite_epsilon_beta_kappa": qmax/(x*N*N) <= h/s and 1/(x*N*N) < 1 and x/N < 1 and x*N*N > 1 and x >= 2,
    "source_target_radius_separation": N/x > 2*rhi and rho > max(F(20), 4*rhi) and 1+rhi < 6,
    "uniform_local_J_and_pin": 9/rho < e < F(1, 2) and Qhi < pin,
    "baseline_log_event": mean < e/2 and 4*variance/e**2 < e/8,
    "two_count_events_and_whole_interval": 4/(e*e*N) < e/8 and 2*h+1/N < 3*e/4,
    "original_target_cap_event": 1+x/(2*N)-1/(2*rhi) < F(19, 20) and target_cap_fail < e/8,
    "new_mu_and_first_profile_event": 0 < mulo < muhi < chi and dlo > 0 and 10/N < dlo/2 and 400/(dlo*dlo*N) < e/8,
    "full_matching_transports": 240/L < e and 120*h < dlo,
    "true_source_face_mean_and_variance": (muhi+2*dhi)/(1-e) < face < chi and s/(1-e) < 10,
    "new_face_retention_at_2pow27": face_gap > 0 and face_fail < 1-Ret,
    "historical_2pow26_failure_preserved": old_face_fail > 1-Ret,
    "new_mass_degree_deficit": (1-e)*masslo > mass > gamma and mass*Ret > gamma and 0 < P0 < 1 and M > 2 and b > 0 and (1-e)*Gamma > gamma,
    "weighted_source_cap_and_tail": nu > 1/N and eta > 0 and K > 4 and cap < e,
    "new_tau_Lsquare_and_tail": ct == F(173, 2000000) > 0 and ct*L > 4 and tau < e and -2-F(4)+F(4)**3/6 > 0,
    "full_H_cross_L_and_tail": sig-dt == F(1, 50000) > 0 and min(sig, dt, sig-dt)*L > 2 and scalar < allowance,
    "price_complete_scalar_allowance": pricehi < price and gamma-price-c == allowance > 0,
    "B_positive_and_same_history_guard": L >= 8 and x <= L and L/4 > gamma and N >= 2 and 8*N > 4 and rho > 6,
}
assert all(checks.values()), checks
values = dict(N12=N, mu_lower=mulo, mu_upper=muhi, dmu_lower=dlo,
              mass_lower_after_error=(1-e)*masslo, paid_Gamma=(1-e)*Gamma,
              baseline_failure=4*variance/e**2, count_failure=4/(e*e*N),
              target_cap_failure=target_cap_fail, first_profile_failure=400/(dlo*dlo*N),
              face_fail_2pow27=face_fail, face_fail_2pow26=old_face_fail,
              weighted_cap_exponent=K, weighted_cap_T8_upper=cap,
              tau_rate=ct, tau_T8_upper=tau, H_cross_rate=sig-dt,
              H_full_remainder_T8_upper=scalar, price_log_upper=pricehi)
receipt = {
    "status": "PASS", "dependency": "R088 accepted generic proof, main acceptance and full review hashes all match",
    "registration_sha256": sha256((HERE/"REGISTRATION.json").read_bytes()).hexdigest(),
    "dependency_acceptance_sha256": sha256((HERE/"DEPENDENCY_ACCEPTANCE.json").read_bytes()).hexdigest(),
    "recipe": reg["recipe"], "checks": checks,
    "exact_values": {k: str(v) for k, v in values.items()},
    "display_only": {k: float(v) for k, v in values.items()},
    "no_producer_program_executed": True, "old_failure_not_overwritten": True,
    "scope": "20 compact independent base/history checks plus analytic tails in REPORT.md; no standalone operator acceptance"
}
(HERE/"verification.json").write_text(json.dumps(receipt, indent=2)+"\n")
print(json.dumps({"status": "PASS", "checks": len(checks),
                  "face27": float(face_fail), "face26_failure": float(old_face_fail),
                  "cap_T8": float(cap), "tau_T8": float(tau), "H_remainder_T8": float(scalar)}))
