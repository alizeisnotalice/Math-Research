#!/usr/bin/env python3
"""New R078 Bernstein gates only; independent rational arithmetic."""
from fractions import Fraction as F
from hashlib import sha256
from math import factorial
from pathlib import Path
import json

HERE = Path(__file__).resolve().parent
SOURCE = HERE.parents[1] / "round073/cantelli_frontier/PROOF.md"
SOURCE_HASH = "50d48b85f02a14f50614f0eea5cedbd7ccd68e9b2f7e1564fe8b41fcd2af9c4b"
assert sha256(SOURCE.read_bytes()).hexdigest() == SOURCE_HASH

x0, chi, p, b, v = F(7392), F(63, 1280), F(183, 1000), F(183, 1000), F(531, 100000)
eB = eN = eT = F(1, 2**19)
eJ = eRad = eS = eD = F(1, 2**20)
m = F(2211, 50000) / (1-eJ)
a = chi-m
u0 = a*x0-p
H0 = u0*u0/(2*(v*x0+b*u0/3))
prefix = sum((F(16)**k/F(factorial(k)) for k in range(33)), F(0))
# log(6/5)=2 atanh(1/11); four terms and exact geometric tail.
z = F(1, 11)
qlo = 2*sum((z**(2*k+1)/F(2*k+1) for k in range(4)), F(0))
qhi = qlo + 2*z**9/(9*(1-z*z))
Ret = F(999999, 1000000)
P0 = 1-eB-F(6, 768)
Ulo = P0*(1-eN)/(1+eJ)
Uhi = (1+eN)/(P0*(1-eJ))
M = 2*(1+eRad)*Uhi
brow = 2*Ret*(1-eRad)*Ulo-eS-eT
Gamma = F(3043, 10000)*brow/M
gamma, oldgamma, price, c = F(599, 2000), F(291, 1000), F(257, 2500), F(393, 2000)
checks = {
    "source_sha256_matches": True,
    "Q0_upper_lt_b": 0 < qlo < qhi < b,
    "positive_a_and_u0": a > 0 and u0 > 0,
    "H0_gt_16": H0 > 16,
    "exp_prefix_gt_inverse_failure": prefix > 10**6,
    "retention_strictly_improves": Ret > F(9717, 10000),
    "new_deficit_gate": 0 < gamma < (1-eD)*Gamma < Gamma < 1,
    "positive_B_minus_gamma_from_old_Lbound": F(8192, 4) > gamma,
    "gamma_difference": gamma-oldgamma == F(17, 2000),
    "same_scalar_remainder": gamma-price-c == F(1, 5000),
}
assert all(checks.values()), checks
values = dict(a=a, u0=u0, H0=H0, exp_prefix=prefix, Q0_lower=qlo,
              Q0_upper=qhi, Ret=Ret, Gamma=Gamma, paid_Gamma=(1-eD)*Gamma,
              gamma=gamma, c=c, full_F_gamma_gain=gamma-oldgamma)
receipt = {
    "status": "PASS", "scope": "10 new rational gates, not the full operator theorem",
    "source_sha256": SOURCE_HASH,
    "checks": checks,
    "exact_values": {k: str(value) for k, value in values.items()},
    "display_only": {k: float(value) for k, value in values.items()},
}
(HERE/"verification.json").write_text(json.dumps(receipt, indent=2)+"\n")
print(json.dumps({"status": "PASS", "new_gates": len(checks),
                  "Gamma": float(Gamma), "gamma": float(gamma), "c": float(c)}))
