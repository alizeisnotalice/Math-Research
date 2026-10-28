#!/usr/bin/env python3
"""Independent compact R080 gates; does not execute producer's 53-gate code."""
from fractions import Fraction as F
from hashlib import sha256
from math import factorial
from pathlib import Path
import argparse
import json

HERE = Path(__file__).resolve().parent
parser = argparse.ArgumentParser()
parser.add_argument("--compute-dir", type=Path, default=HERE.parent/"compute")
args = parser.parse_args()
reg = json.loads((HERE/"REGISTRATION.json").read_text())
for name, digest in reg["producer_locks"].items():
    assert sha256((args.compute_dir/name).read_bytes()).hexdigest() == digest
source = HERE.parents[1]/"round073/cantelli_frontier/PROOF.md"
assert sha256(source.read_bytes()).hexdigest() == reg["frozen_R073_proof_sha256"]
producer = json.loads((args.compute_dir/"case_988_1000.json").read_text())

def log12(y):
    assert y > 1
    z = (y-1)/(y+1)
    lo = 2*sum((z**(2*j+1)/F(2*j+1) for j in range(12)), F(0))
    return lo, lo+2*z**25/(25*(1-z*z))

def S32(y):
    assert y > 0
    return sum((y**j/F(factorial(j)) for j in range(33)), F(0))

P = {k: F(v) for k, v in reg["tuple"].items()}
L, sig, chi, qlo, qhi, D = (P[k] for k in ("L0", "sigma", "chi", "qmin", "qmax", "hden"))
x = (1-sig)*L
N = L**32/F(factorial(32))
h = 1/(D*L)
rlo, rhi = 1+P["trim"], F(3, 2)*qhi-P["trim"]
eB = eN = eT = F(1, 2**19)
eJ = eR = eS = eD = F(1, 2**20)
ea = F(1, 2**16)
P0 = 1-eB-6/D
Ulo = P0*(1-eN)/(1+eJ)
Uhi = (1+eN)/(P0*(1-eJ))
M = 2*(1+eR)*Uhi
brow = 2*P["Ret"]*(1-eR)*Ulo-eS-eT
Gamma = P["mL"]*brow/M
mlog, _ = log12(rhi/rlo)
Qlo, Qhi = log12(F(6, 5))
mu1hi = F(1, 2)-F(5, 2)*Qlo
mu2lo = 1-5*Qhi-F(5, 2)*Qhi**2
mu2hi = 1-5*Qlo-F(5, 2)*Qlo**2
H2lo = (1-1/N)*(mu2lo-(mu2hi+F(1, 25))*h)
H2hi = mu2hi+1/(250*N*qlo)
ML = F(3, 4)*L**2/(qlo*qlo*N)+F(1, 2)*L/(qlo*N)
VL = F(3, 2)*L**2/(qlo*qlo*N)
pin, m, v = F(183, 1000), F(2211, 50000)/(1-eJ), F(531, 100000)
a = chi-m
u = a*x-pin
H = u*u/(2*(v*x+pin*u/3))
eta = (1-qhi)/(2*qhi)
nu = F(2, 3)-F(5, 3)*eN
K = (nu-1/N)*eta*x/(1+eJ)
cT = sig/2-chi*(1-sig)
cap = M/S32(K)
tau = rhi*D*L*64/S32(cT*L)
root = 1/N+1/S32(sig*L)+128/S32(sig*L/2)
_, pricehi = log12(1/(1-sig))

checks = {
    "producer_identity_not_reexecuted": producer["status"] == "PASS" and len(producer["gates"]) == 53 and all(producer["gates"].values()),
    "exact_selected_tuple": x == F(91136, 25) and h == F(1, 6291456) and rhi == F(741, 500)-F(1, 2**19),
    "active_band_and_radius": qlo == F(2, 3) and qlo < qhi < 1 and 1 < rlo < rhi < F(3, 2) and 2*h < P["trim"] and h < F(1, 1000) and N*qlo/x > 4*rhi and N/x > rhi,
    "epsilon_beta_kappa_guard": qhi/(x*N*N) <= h/F(3, 2) and 1/(x*N*N) < F(1, 2) and x/N < 1 and x*N*N > 1 and 2 <= x <= L,
    "baseline_log_probability_budget": ML <= eB/2 and 4*VL/eB**2 <= ea/8,
    "both_count_and_transport": 4/(eN**2*N) <= ea/8 and 2*h+1/N <= 3*eN/4 and x/(3*qlo*N) <= eJ,
    "target_cap_probability_budget": 1+x/(2*N)-1/(2*rhi) <= F(3, 4) and N >= 75 and F(3, 4)**74 <= ea/8,
    "first_profile_probability_and_transport": mu1hi < F(221, 5000) and N >= 9000 and F(324000000)/N <= ea/8 and h <= F(1, 100000),
    "second_profile_mean_and_probability": 0 < mu2lo < mu2hi and H2lo > F(528, 100000) and H2hi < F(5295, 1000000) and F(5760000)/N <= ea/8 and h/5 < F(1, 1000000),
    "local_H2_and_true_variance_guard": h <= F(1, 2**20) and N >= 2000 and x/N <= F(1, 1000) and F(5301, 1000000)/(1-eJ) < v and Qhi < pin,
    "new_Bernstein_prefix": a > 0 and u > 0 and S32(H) > 10**6,
    "new_mass_and_degree": 0 < P0 < 1 and M > 2 and brow > 0 and (1-ea)*mlog > P["mL"] and 1/(1-(h/rlo)**2) <= 1+eR,
    "new_weighted_cap_prefix": N > 1/nu and K > 4 and cap < eS,
    "new_tau_prefix": cT == F(961, 100000) and cT*L > 4 and tau < eT and 64**2 == L,
    "new_gamma_budget_and_positive_tau": 0 < P["gamma"] < (1-eD)*Gamma < Gamma < 1 and L/4 > P["gamma"] and (L/2)**2/2 >= L,
    "new_price_and_full_root_remainder": pricehi < P["price"] and sig*L/2 > 4 and root < P["scalar_slack"] and P["gamma"]-P["price"]-P["c"] == P["scalar_slack"],
    "full_cost_vs_Phi69_on_its_domain": P["c"] > F(317, 1000)-F(83, 1000)+F(4, 2**30),
    "full_cost_vs_Phi78_on_its_domain": P["c"] > F(599, 2000)-F(25, 256)+F(4, 8192),
}
assert all(checks.values()), checks
values = dict(Nstar=N, x0=x, h0=h, P0=P0, Ulo=Ulo, Uhi=Uhi,
              M=M, brow=brow, Gamma=Gamma, paid_Gamma=(1-eD)*Gamma,
              mL=P["mL"], gamma=P["gamma"], c=P["c"],
              Bernstein_a=a, Bernstein_gap=u, Bernstein_exponent=H,
              weighted_cap_exponent=K, weighted_cap_upper=cap,
              tau_upper=tau, root_remainder=root, price_upper=pricehi,
              Phi69_comparison_margin=P["c"]-F(234, 1000)-F(4, 2**30),
              Phi78_comparison_margin=P["c"]-(F(599, 2000)-F(25, 256))-F(4, 8192))
receipt = {
    "status": "PASS", "recipe": reg["independent_recipe"],
    "registration_sha256": sha256((HERE/"REGISTRATION.json").read_bytes()).hexdigest(),
    "producer_compute_dir": str(args.compute_dir), "checks": checks,
    "exact_values": {k: str(val) for k, val in values.items()},
    "display_only": {k: float(val) for k, val in values.items()},
    "scope": "18 compact independent gates; mathematical proof and full-chain independent acceptance remain separate"
}
(HERE/"verification.json").write_text(json.dumps(receipt, indent=2)+"\n")
print(json.dumps({"status": "PASS", "new_compact_gates": len(checks),
                  "gamma": str(P["gamma"]), "c": str(P["c"]),
                  "cap_upper": float(cap), "tau_upper": float(tau)}))
